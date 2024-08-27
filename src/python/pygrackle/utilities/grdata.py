#!/usr/bin/env python3

# A tool for managing grackle data files. This should be usable as
# -> a standalone command line tool (when pygrackle isn't installed)
# -> as a part of pygrackle.


# Notes on the file registry
# --------------------------
# The file registry refers to a small file that associates a filename with a checksum.
#
# In the long term, we plan to support 3 cases involving hashes checksums:
#
# 1. have the C layer of Grackle provide the option to directly read in files from the
#    data-store without being provided a full path
#     -> this logic is already mostly implemented
#     -> In this case, we will also need to have access to the file checksum so we can
#        validate that the correct file is being accessed. (This is mostly to ensure we
#        don't break other people's results because we make a mistake). The checksum
#        validation will be performed with a tool like picohash
#     -> When we do this, we will directly embed the information encoded in the file
#        registry inside of scikit-build-core (we picked the file format to ensure that
#        the information is easy to embed in a C file)
#
# 2. Continue supporting the functionality (and cli) implemented by this file within
#    pygrackle
#
# 3. Support running using this script as a standalone command-line program.
#    -> See the end of this file for what that entails


# With that said, the pooch package seems like a great tool for helping us do
# this (it does a bunch of the things we would probably need to do anyways).
# If/when we want to distribute this file as a standalone application alongside
# the Grackle library, we have 2 options:
#  1. We could remove pygrackle as a dependency.
#  2. We could package this file, pooch and pooch's dependencies as a
#     `zipapp <https://docs.python.org/3/library/zipapp.html#module-zipapp>`_.
#     This is only possible as long as:
#       - pooch doesn't add many more required dependencies (right now it only
#         has 2 dependencies)
#       - pooch and its dependencies remain written in pure python
#       - pooch and its dependencies maintain compatability with a wide range of
#         python versions


_DESCRIPTION = """\
This is a command line tool for downloading and managing data files to be used
with Grackle.



"""

import argparse
from contextlib import contextmanager, ExitStack
import filecmp
import hashlib
import io
import os
import re
import shutil
import stat
import sys
import traceback
import warnings

import pooch  # external import

# pooch will assume that any unlabeled checksum was computed with the following algorithm
_POOCH_DEFAULT_CKSUM_KIND = "sha256"

if (sys.version_info.major, sys.version_info.minor) < (3, 3):
    raise RuntimeError("python 3.3 or newer is required")


_UNSPECIFIED = object()
_OBJECT_STORE_SUBDIR = "object-store"


class ToolConfig:
    """Tracks basic information about this tool"""

    def __init__(self, *, grackle_version, protocol_version="1", checksum_kind="sha1"):
        self.grackle_version = grackle_version
        self.protocol_version = protocol_version
        self.checksum_kind = checksum_kind


class DataStoreConfig:
    """Track basic configuration information

    In principle, this information is intended to be a little more
    flexible and might not be known as early as ToolConfig.
    """

    def __init__(
        self,
        *,
        data_dir,
        store_location,
        data_repository_url,
        contemporaneous_git_hash,
        checksum_kind,
        file_registry_file,
    ):
        # where data is actually stored
        self.data_dir = data_dir
        self.store_location = store_location

        # properties for tracking files
        self.data_repository_url = data_repository_url
        self.contemporaneous_git_hash = contemporaneous_git_hash

        # specifies the kinds of checksums listed in the registry
        self.checksum_kind = checksum_kind
        # the following specifies the file containing the file registry
        self.file_registry_file = file_registry_file


def _get_platform_data_dir(appname="grackle"):
    """Returns a string specifying the default data directory

    All of these choices are inspired by the API description of the platformdirs python
    package
        * we only looked at online documentation:
          https://platformdirs.readthedocs.io/en/latest/
        * we have NOT read any source code
    """
    if sys.platform.startswith("win32"):
        raise RuntimeError()
    elif sys.platform.startswith("darwin"):
        return os.path.expanduser(f"~/Library/Application Support/{appname}")
    else:  # assume linux/unix
        # https://specifications.freedesktop.org/basedir-spec/latest/
        dflt = "~/.local/share"
        env_str = os.getenv("XDG_DATA_HOME", default=dflt)
        if env_str[:1] not in ["~", "/"]:
            # this is what the specification tells us to do
            warnings.warn(
                "ignoring XDG_DATA_HOME because it doesn't hold an " "absolute path"
            )
            env_str = dflt

        # now actually infer the absolute path
        if (env_str[0] == "~") and (not env_str.startswith("~/")):
            if env_str[:2] != "~/":  # for parity with C-version of this function
                raise RuntimeError(
                    "can't expand can't expand env-variable, XDG_DATA_HOME when "
                    "it starts with `~user/` or just contains `~`"
                )
            return os.path.expanduser(f"{env_str}/{appname}")
        else:
            return f"{env_str}/{appname}"


def _get_data_dir():
    manual_choice = os.getenv("GRACKLE_DATA_DIR", default=None)
    if (manual_choice is None) or (len(manual_choice) == 0):
        return _get_platform_data_dir()
    elif (manual_choice[0] != "~") and (not os.path.isabs(manual_choice)):
        raise RuntimeError("GRACKLE_DATA_DIR must specify an absolute path")
    elif manual_choice[0] == "~":
        if not env_str[:2] != "~/":  # for parity with C-version of this function
            raise RuntimeError(
                "can't expand can't expand env-variable, GRACKLE_DATA_DIR when "
                "it starts with `~user/` or just contains `~`"
            )
        return os.path.expanduser(manual_choice)
    else:
        return manual_choice


def _parse_file_registry(f):
    """Read the file registry, as a dict from a text file

    Parameters
    ----------
    f : file or str or bytes or ``os.PathLike``
        Contains the data to be read in

    Notes
    -----
    We describe the format below. This format was choosen so that the
    contents could be injected into a C to be used as a literal.

      * empty lines and lines that start with ``//`` are ignored

      * all other lines should look like ``{"<file-name>", "<hash>"}``
        and there is allowed to be a trailing comma
    """

    if (sys.version_info.major, sys.version_info.minor) < (3, 6):
        if not isinstance(f, io.IOBase):
            path = f
        else:
            path = None
    else:
        try:
            path = os.fspath(f)
        except TypeError:
            path = None

    with ExitStack() as stack:
        if path is None:
            file = f
        else:
            file = stack.enter_context(open(path, "r"))

        file_registry = {}
        for i, line in enumerate(file):  # iterater over lines
            if (len(line) == 0) or line.isspace() or line.startswith("//"):
                continue
            m = re.match(
                r'^\s*{\s*"(?P<fname>[^"]+)"\s*,\s*"(?P<cksum>[^"]+)"\s*},?\s*', line
            )
            if m is None:
                raise RuntimeError(
                    f"Something went wrong with parsing line {i+1} of {f}:\n "
                    f"  `{line}`"
                )
            file_registry[m["fname"]] = m["cksum"]
    return file_registry


class GenericToolError(RuntimeError):
    pass


class LockFileExistsError(FileExistsError):
    pass


@contextmanager
def lock_dir(lock_file_path):
    """
    Contextmanager that creates a "lock file." The context-manager will delete
    the file when we finish. If the lock already exists, the program will abort
    with an explanatory error (this ensures that only 1 copy of the program will
    try to run at a time).
    """
    try:
        f = open(lock_file_path, "x")
        f.close()
    except FileExistsError as err:
        raise LockFileExistsError(
            err.errno,
            err.strerror,
            err.filename,
            getattr(err, "winerror", None),
            err.filename2,
        ) from None

    try:
        yield None
    finally:
        os.remove(lock_file_path)


@contextmanager
def standard_lockfile(data_config):
    lock_file_path = os.path.join(data_config.data_dir, "lockfile")
    with lock_dir(lock_file_path):
        yield None


def calc_checksum(fname, alg_name, *, chunksize=4096):
    """Calculate the checksum for a given fname"""
    # construct the object to track intermediate state of the checksum
    # calculation as we stream through the data
    hash_obj = hashlib.new(alg_name)
    with open(fname, "rb") as f:
        buffer = bytearray(chunksize)
        while True:
            nbytes = f.readinto(buffer)
            if nbytes == chunksize:
                hash_obj.update(buffer)
            elif nbytes:  # equivalent to: (nbytes is not None) and (nbytes > 0)
                hash_obj.update(buffer[:nbytes])
            else:
                break
    return hash_obj.hexdigest()


def matches_checksum(fname, alg_name, checksum):
    return checksum == calc_checksum(fname, alg_name)


def _create_retriever(destination_path, data_config):
    """
    create pooch object responding for fetching data files

    Notes
    -----
    If we ever move away from pooch (e.g. to make this functionality easier to use in a
    portable standalone script), we need to ensure that the our new approach implements
    a similar procedure that they adopt where:
        1. any downloaded file is first put in a temporary location
        2. and then, only after we verify that the checksum is correct, we move the
           file to the downloads directory.
    """
    repo_url = data_config.data_repository_url
    repo_version = data_config.contemporaneous_git_hash
    prefix = f"{data_config.checksum_kind}"

    file_registry = _parse_file_registry(data_config.file_registry_file)

    # if we move away from pooch, (
    return pooch.create(
        path=destination_path,
        base_url=f"{repo_url}/raw/{repo_version}/input/",
        registry=dict((k, f"{prefix}:{v}") for k, v in file_registry.items()),
    )


def _pretty_log(arg):
    """indent messages so it's clear when multiline messages are a single thought"""
    lines = arg.splitlines()
    if len(lines):
        print("\n".join([f"-- {lines[0]}"] + [f"   {e}" for e in lines[1:]]))
    else:
        print("")


def _ensure_exists(path, content_description):
    if not os.path.isdir(path):
        if len(content_description) > 0:
            _pretty_log(f"creating directory {content_description}\n-> {path}")
        os.mkdir(path)


# to be used with os.chmod to set permissions to prevent mutations of files (you can
# always delete it if you own it)
_IMMUTABLE_MODE = stat.S_IREAD | stat.S_IRGRP | stat.S_IROTH


class _HardlinkStrat:
    """
    Acts as a "namespace" for functions related to our deduplication strategy
    that uses Hardlinks
    """

    @staticmethod
    def are_linked(fname, fname2):
        """return whether ``fname`` & ``fname2`` specify paths that are hardlinks"""
        try:
            statinfo1 = os.stat(fname, follow_symlinks=False)
            statinfo2 = os.stat(fname2, follow_symlinks=False)
        except FileNotFoundError:
            return False
        return statinfo1.st_ino == statinfo2.st_ino

    @staticmethod
    def remove_if_norefs(fname):
        """
        Removes the specified file if there are no other references to it.

        Parameters
        ----------
        fname : str
            Path to the file that we are operating on

        Returns
        -------
        bool
            Indicates if any file was removed
        """
        statinfo = os.stat(fname, follow_symlinks=False)

        # statinfo.st_nlink == 1 means that the only hardlink is the hardlink of
        # associated with fname
        # -> it should be possible for ``os.stat(fname).st_nlink`` to return ``0``
        if statinfo.st_nlink == 1:
            os.remove(fname)
            return True
        return False

    @staticmethod
    def deduplicate(full_fname, shared_fname):
        """
        Perform logic to ensure that ``full_fname`` and ``shared_fname``
        are both paths that refer to the same hardlink.

        This handles 3 main cases:

          1. ``full_fname`` and ``shared_fname`` are already hardlinked.

             * Nothing is done.

          2. ``full_fname`` exists and ``shared_fname`` doesn't.

             * A hardlink will be created at ``shared_fname`` that refers
               to ``full_fname``.

          3. ``full_fname`` and ``shared_fname`` specify existing distinct
             copies of the same existing file.

             * in this case, ``full_fname`` is deleted and then replaced
               with a hardlink that refers to ``shared_fname``.

        Parameters
        ----------
        full_fname : str
            Specifies an existing file-path that already exists
        shared_fname : str
            Specifies a file-path that may or may not exist. If it does
            already exist, it will be preserved (in case it is already in
            use for deduplicating other existing files)
        """
        if not os.path.isfile(full_fname):
            raise FileNotFoundError(full_fname)
        elif _HardlinkStrat.are_linked(full_fname, shared_fname):
            pass  # do nothing!
        elif os.path.isfile(shared_fname):
            if not filecmp.cmp(full_fname, shared_fname, shallow=False):
                raise ValueError(
                    f"`{full_fname}` and `{shared_fname}` specify files that aren't "
                    "perfect copies"
                )
            os.remove(full_fname)
            os.link(shared_fname, full_fname)
        else:
            os.link(full_fname, shared_fname)


def _fetch_files(retriever, cksum_kind, object_dir):
    """
    Does the heavy lifting of fetching files.

    Parameters
    ----------
    retriever : ``pooch.Pooch``
        pooch object that manages downloads of the specified files
    cksum_kind : str
        The primary checksum algorithm that the tool is globally configured to use.
    object_dir : str
        Path to the object directory. This is the name where checksum names are used as
        filenames. (This is the mechanism used to aid deduplication)
    """
    try:
        import tqdm

        progressbar = True
    except:
        progressbar = False

    num_download_attempts = 0

    # NOTE: the docstring for ``pooch.create`` makes it clear that ``retriever.registry``
    #       is a part of the public API
    for fname, full_checksum_str in retriever.registry.items():
        # extract the checksum_kind and string that are stored in the registry
        # (we are being a little more careful here than necessary, but if this ever
        # becomes library-code, it will pay off)
        if ":" in full_checksum_str:
            cur_cksum_kind, checksum = full_checksum_str.split(":")
        else:
            cur_cksum_kind, checksum = _POOCH_DEFAULT_CKSUM_KIND, full_checksum_str

        if cur_cksum_kind != cksum_kind:
            raise ValueError(
                "Currently, we only support downloading from file registries where the "
                f"checksum algorithm matches the globally used algorithm, {cksum_kind}. "
                f"The checksum algorithm associated with {fname} is {cur_cksum_kind}"
            )

        # name associated with current file in the current grackle version
        full_fname = os.path.join(retriever.path, fname)

        # if the file already exists we are done
        if os.path.exists(full_fname):
            if not matches_checksum(full_fname, cksum_kind, checksum):
                raise RuntimeError(
                    f"{full_fname} already exists but has the wrong hash"
                )
            continue

        num_download_attempts += 1

        # download the file (pooch will log a detailed message
        retriever.fetch(fname, progressbar=progressbar)
        os.chmod(full_fname, _IMMUTABLE_MODE)

        # now deduplicate
        checksum_fname = os.path.join(object_dir, checksum)

        try:
            _HardlinkStrat.deduplicate(full_fname, checksum_fname)

            # not strictly necessary, but doing this for safety reasons
            os.chmod(checksum_fname, _IMMUTABLE_MODE)

        except Exception as err:
            # remove full_fname since we don't want users to use it before dealing
            # with the larger issue. We also want to make the errors reproducible
            os.remove(full_fname)
            if not (isinstance(err, ValueError) and os.path.is_file(checksum_fname)):
                raise err

            # this should only happens when full_fname and checksum_fname both exist,
            # but aren't perfect matches of each other. We try to provide a more
            # informative error message
            if not matches_checksum(
                checksum_fname, data_config.checksum_kind, checksum
            ):
                raise GenericToolError(f"""\
A file (used for deduplication) that already existed on disk
   `{checksum_fname}`
which is probably a version of `{fname}`,
doesn't have the appropriate {data_config.checksum_kind} checksum.
-> expected: {calc_checksum(checksum_fname, data_config.checksum_kind)}
-> actual: {checksum}
-> This implies that the data was corrupted and it needs to be dealt with.
   To avoid confusion we have deleted the newly downloaded version of
   `{fname}`
-> The safest bet is probably to delete the data directory""")
            else:
                raise GenericToolError(f"""\
Something bizare (& extremely unlikely) happened:
-> a previous invocation of this tool appears to have installed a data file
   with the same checksum as {fname}, but has different contents.
-> we adopt a similar system to git and the odds for this to organically
   happen for a small collection of files is truly astronomical!
-> this is probably a sign that something went wrong. We deleted the newly
   downloaded version of the file""")

    if num_download_attempts == 0:
        _pretty_log("no files needed to be loaded")


def get_version_dir(tool_config, data_config):
    return os.path.join(data_config.store_location, tool_config.grackle_version)


def fetch_command(args, tool_config, data_config):
    # the data_dir is a directory that contains:
    # -> the data-store directory for data managed by the current protocol version
    # -> (possibly) data-store directories for data managed by other protocol version
    # -> (possibly) a directory called `user-data/` where users can put custom data
    _ensure_exists(data_config.data_dir, "that will hold all Grackle data")

    with standard_lockfile(data_config):
        # even though it isn't used for anything right now, make the directory that is
        # reserved for user content
        _ensure_exists(
            os.path.join(data_config.data_dir, "user-data"),
            "reserved for user-defined data",
        )

        # do a little more setup!
        _ensure_exists(data_config.store_location, "that will hold the data-store")

        # ensure version_dir and object_dir both exist. They respectively store
        # filenames that are (hard) linked to the data-file.
        # -> version_dir uses the names known by the associated grackle-version
        # -> object_dir uses the checksum as a filenames
        object_dir = os.path.join(data_config.store_location, _OBJECT_STORE_SUBDIR)

        _ensure_exists(object_dir, "")
        version_dir = get_version_dir(tool_config, data_config)
        _ensure_exists(version_dir, "that holds data for current Grackle version")

        # create the object that is used to actually loads the data
        retriever = _create_retriever(
            destination_path=version_dir, data_config=data_config
        )

        _fetch_files(retriever, data_config.checksum_kind, object_dir)


def direntry_iter(path, *, ftype="file", mismatch="skip", ignore=None):
    """
    Iterate over the contents of a single directory with focus on a
    particular file type assumption that all.

    Parameters
    ----------
    path : str
        path to the directory
    ftype : {None, 'file', 'dir'}
        When not ``None``, the iterator only produces entries for the
        specified file-type
    mismatch : {'skip', 'lazy_err', 'eager_err'}
        Specifies the action to take when this generator encounters an
        entry in ``path`` that doesn't have the specified type.
          * ``'skip'`` means to simply skip the entry
          * ``'lazy_err'`` means that we raise an error
          * ``'eager_err'`` means that we check for any mismatches and
            raise an error if any mismatch is encountered and afterwards,
            we start yielding elements
    ignore : container of str, optional
        Optional container of strings that are ignored

    Yields
    ------
    pair : tuple of two str
        The first element is the entry's base filename and the second is the full path
    """
    if ftype is None:
        has_ftype = lambda full_path: True
    if ftype == "dir":
        has_ftype = os.path.isdir
    elif ftype == "file":
        has_ftype = os.path.isfile
    else:
        raise ValueError("ftype must be None, 'file' or 'dir'")

    if ignore is None:
        ignore = []
    elif isinstance(ignore, str):
        raise TypeError("ignore can't be a string")

    it = map(
        lambda e: (e, os.path.join(path, e)),
        filter(lambda e: e not in ignore, os.listdir(path)),
    )
    if mismatch == "eager_err":
        for pair in direntry_iter(path, ftype=ftype, mismatch="lazy_err"):
            pass
        yield from it
    elif mismatch in ["lazy_err", "skip"]:
        for pair in it:
            if has_ftype(pair[1]):
                yield pair
            elif mismatch == "lazy_err":
                raise RuntimeError(f"{pair[1]} isn't a {ftype}")
    else:
        raise ValueError("mismatch must be 'eager_err', 'lazy_err' or 'skip'")


def rm_command(args, tool_config, data_config):
    """Logic for removing files"""
    if args.vdata is _UNSPECIFIED:
        # this means that we are removing the whole data store
        if not args.data_store:
            raise RuntimeError("SOMETHING WENT HORRIBLY, HORRIBLY WRONG")

        _descr = os.path.basename(data_config.store_location)
        target_path = data_config.store_location
        operation_description = (
            f"deleting ALL files in the data-store associated with this tool, {_descr}"
        )
        if not os.path.isdir(target_path):
            raise GenericToolError(
                "intended to recursively delete all contents of the associated data-store. "
                "But no such directory can be found."
            )

        fn = shutil.rmtree

    else:
        if args.vdata is None:
            target = tool_config.grackle_version
            _descr = f"associated with this tool (`{tool_config.grackle_version}`)"
        else:
            target = args.vdata
            _descr = f"`{target}`"
        target_path = os.path.join(data_config.store_location, target)
        operation_description = (
            f"deleting all data file references for the grackle-version {_descr}. "
            "Any files for which the reference-count drops to zero will also be removed."
        )

        if not os.path.isdir(target_path):
            raise GenericToolError(
                "intended to delete all data-file references for the grackle-version "
                f"{_descr}, but no such data is tracked in the data-store."
            )

        def fn(path):
            object_dir = os.path.join(data_config.store_location, _OBJECT_STORE_SUBDIR)
            if not os.path.isdir(object_dir):
                raise RuntimeError(
                    "SOMETHING IS HORRIBLY WRONG!!! THE {object_dir} IS MISSING"
                )

            # we throw an err if this directory contains some unexpected stuff
            it = direntry_iter(path, ftype="file", mismatch="eager_err")
            for name, full_path in it:
                # get path to corresponding hardlinked file in _OBJECT_STORE_SUBDIR
                checksum = calc_checksum(full_path, alg_name=tool_config.checksum_kind)
                checksum_fname = os.path.join(object_dir, checksum)
                checksum_fname_exists = os.path.isfile(checksum_fname)

                if not checksum_fname_exists:
                    warnings.warn(
                        "Something weird has happened. There is no deduplication file "
                        f"associated with {full_path}"
                    )
                os.remove(full_path)
                if checksum_fname_exists:
                    _HardlinkStrat.remove_if_norefs(checksum_fname)
            os.rmdir(path)

    target_exists = os.path.isdir(target_path)

    if not args.force:
        _pretty_log(
            f"{operation_description}\n"
            "-> essentially, we are recursively removing\n"
            f"     `{target_path}`\n"
            "-> to actually perform this command, pass the --force flag"
        )
    else:
        fn(target_path)


def lsversions_command(args, tool_config, data_config):
    if not os.path.exists(data_config.store_location):
        print("there is no data")
    with standard_lockfile(data_config):
        it = direntry_iter(
            data_config.store_location,
            ftype="dir",
            mismatch="lazy_err",
            ignore=[_OBJECT_STORE_SUBDIR],
        )
        print(*sorted(pair[0] for pair in it), sep="\n")


def getpath_command(args, tool_config, data_config):
    print(data_config.data_dir)


def calcreg_command(args, tool_config, data_config):
    # print the properly file registry information (in the proper format that can be
    # used to configure newer versions of Grackle

    # we use listdir since we are targetting 3.3, but we set things up so that we could
    # use os.scandir
    try:
        it = direntry_iter(args.path, ftype="file", mismatch="eager_err")
    except FileNotFoundError:
        raise ValueError(f"{path!r} doesn't specify a directory or file")
    except NotADirectoryError:
        it = [(os.path.basename(args.path), args.path)]

    pairs = [(name, calc_checksum(path, args.hash_name)) for name, path in it]

    with ExitStack() as stack:
        if args.output is None:
            file = sys.stdout
        else:
            file = stack.enter_context(open(args.output, "w"))

        if file is not None:
            file.write(f"""\
// This is a file registry generated by the grackle data management tool
// To overwrite this file with an updated copy (assuming that pygrackle is
// installed), you might invoke:
//    python -m pygrackle --hash_name {args.hash_name} --output <outpath> <dir>
// in this sample command, you would substitute:
//    -> ``<outpath>`` with a path to the output file
//    -> ``<dir>`` with a path to the directory containing all files that are
//       to be included in the registry
""")
        print(*[f'{{"{p[0]}", "{p[1]}"}}' for p in sorted(pairs)], sep=",\n", file=file)


def _add_version(parser, version_flag, version_name, value):
    """add argument to parser to show a version and exit (similar to --help)"""

    class _Action(argparse.Action):
        def __call__(self, *args, **kwargs):
            print(value)
            sys.exit(0)

    parser.add_argument(
        version_flag,
        metavar="",
        action=_Action,
        nargs=0,
        help=f"show associated {version_name} and exit",
    )


def build_parser(tool_config, prog_name):
    parser = argparse.ArgumentParser(
        prog=prog_name,
        description=_DESCRIPTION,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    _add_version(
        parser, "--version-grackle", "Grackle version", tool_config.grackle_version
    )
    _add_version(
        parser,
        "--version-protocol",
        "data-store protocol version",
        tool_config.protocol_version,
    )

    subparsers = parser.add_subparsers(required=True)

    parser_fetch = subparsers.add_parser(
        "fetch",
        help=(
            "fetch data files if we don't already have the data for the "
            "associated version of grackle"
        ),
    )
    parser_fetch.set_defaults(func=fetch_command)

    parser_ls = subparsers.add_parser("ls-versions", help="list the versions")
    parser_ls.set_defaults(func=lsversions_command)

    parser_rm = subparsers.add_parser(
        "rm", help="remove data associated with a given version"
    )
    parser_rm.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="This option must be present to actually remove things",
    )
    rm_spec_grp = parser_rm.add_argument_group(
        title="Target", description="specifies the target that will be removed"
    ).add_mutually_exclusive_group(required=True)
    rm_spec_grp.add_argument(
        "--data-store", action="store_true", help="remove the full data-store"
    )
    rm_spec_grp.add_argument(
        "--vdata",
        default=_UNSPECIFIED,
        nargs="?",
        help="remove all data associated with the contemporaneous grackle version",
    )
    parser_rm.set_defaults(func=rm_command)

    parser_getpath = subparsers.add_parser(
        "getpath", help="get filesystem location where all the data is stored"
    )
    parser_getpath.set_defaults(func=getpath_command)

    parser_calcregistry = subparsers.add_parser(
        "calcreg",
        help=(
            "prints the file registry (file hash pairs) for a given directory. This "
            "computed registry can be used to configure future versions of Grackle."
        ),
    )
    parser_calcregistry.add_argument(
        "-o",
        "--output",
        metavar="FILE",
        help=(
            "Write the output to a file instead of stdout. The file will include extra "
            "metadata (as comments)."
        ),
    )
    parser_calcregistry.add_argument(
        "--hash-name",
        required=True,
        metavar="HASH",
        choices=hashlib.algorithms_guaranteed,
        help=(
            "the kind of checksum to compute. Must be one of: "
            f"{ ', '.join(sorted(hashlib.algorithms_guaranteed))}"
        ),
    )
    parser_calcregistry.add_argument(
        "path", help="path to the directory containing the files in the registry"
    )
    parser_calcregistry.set_defaults(func=calcreg_command)

    return parser


def main(tool_config, data_config, prog_name):
    parser = build_parser(tool_config, prog_name)
    args = parser.parse_args()

    try:
        args.func(args, tool_config=tool_config, data_config=data_config)
    except SystemExit:
        pass  # this shouldn't come up!
    except LockFileExistsError as err:
        lock_file_path = err.filename
        print(
            f"""\
ERROR: The `{lock_file_path}` lock-file already exists.
-> This probably means that another copy of this tool is currently running.
-> If you are absolutely sure that's not the case, that probably means that a copy
   of this tool previously crashed""",
            file=sys.stderr,
        )
        sys.exit(78)  # https://www.man7.org/linux/man-pages/man3/sysexits.h.3head.html
    except GenericToolError as err:
        print(f"ERROR: {err.args[0]}")
        sys.exit(70)  # https://www.man7.org/linux/man-pages/man3/sysexits.h.3head.html
    except:
        print(f"Unexpected error:", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        sys.exit(70)  # https://www.man7.org/linux/man-pages/man3/sysexits.h.3head.html
    else:
        sys.exit(0)


def _default_data_config(tool_config, file_registry_file):
    """Provides default data configuration"""
    _REPO_URL = "https://github.com/grackle-project/grackle_data_files/"

    # this is hash that holds the versions of the datafiles from the time when this
    # version of the file was shipped
    _CONTEMPORANEOUS_COMMIT_HASH = "9a63dbefeb1410483df0071eefcbff666f40816d"

    # FILE_REGISTRY is in a format that could be injected into a C file as a literal
    data_dir = _get_data_dir()
    protocol_version = tool_config.protocol_version
    return DataStoreConfig(
        data_dir=data_dir,
        store_location=os.path.join(data_dir, f"data-store-v{protocol_version}"),
        data_repository_url=_REPO_URL,
        contemporaneous_git_hash=_CONTEMPORANEOUS_COMMIT_HASH,
        checksum_kind=tool_config.checksum_kind,
        file_registry_file=file_registry_file,
    )


def make_config_objects(grackle_version, file_registry_file):
    """Construct the pair of configuration objects used for running the calculation

    Parameters
    ----------
    grackle_version : str
        the version of grackle (NOT pygrackle)
    file_registry_file : file or str or bytes or ``os.PathLike``
        Contains the file registry
    """
    tool_config = ToolConfig(grackle_version=grackle_version)
    data_config = _default_data_config(tool_config, file_registry_file)
    return tool_config, data_config


# to support installing this as a standalone script, we will need to introduce the
# following procedure to the build-system:
#    - treat this file as a template-file and configure it with CMake's
#      ``configure_file`` command (or invoke ``configure_file.py`` under the classic
#      build system) in order to substitute the names enclosed by the @ symbols
#    - if we are still using pooch (or some other external package) we'll need to
#      introduce logic to convert this into a zipapp (this is a special zip file that
#      contains all dependencies that the python interpretter knows how to execute)
#    - make resulting file executable (and maybe drop the .py suffix)
#    - install it into the bin directory alongside the grackle libraries

if __name__ == "__main__":
    _GRACKLE_VERSION = "@GRACKLE_VERSION@"
    _FILE_REGISTRY_CONTENTS = """\
@FILE_REGISTRY_CONTENTS@
"""

    def _check_substitution_problems(var_name, var_value):
        if (
            (var_name in var_value)
            or ("@" in var_value)
            or (len(var_value) == 0)
            or (var_value.isspace())
        ):
            raise RuntimeError(
                "something went wrong when the build-system was configuring the "
                f"{var_name} variable"
            )

    _check_substitution_problems("GRACKLE_VERSION", _GRACKLE_VERSION)
    _check_substitution_problems("FILE_REGISTRY_CONTENTS", _FILE_REGISTRY_CONTENTS)

    _CONFIG_PAIR = make_config_objects(
        grackle_version=_GRACKLE_VERSION,
        file_registry_file=io.StringIO(_FILE_REGISTRY_CONTENTS),
    )
    main(*_CONFIG_PAIR, prog_name="grdata")
