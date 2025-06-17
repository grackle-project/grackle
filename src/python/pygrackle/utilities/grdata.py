#!/usr/bin/env python3
"""A tool for managing grackle data files usable as part of pygrackle and as a
standalone command line tool (when pygrackle IS NOT installed)

Just below, we provide a detailed description that will also be used for documentation
on the website AND as queryable documentation via the ``help`` subcommand. We enclose
some text in triple braces for 2 purposes:
- it is designed to be anchors used by sphinx to include the documentation.
- to help format output of the ``help`` subcommand. Anchors of the format
  ``[[[BEGIN-SECTION:<name>]]]`` will be replaced with a section title ``<name>``,
  and all other anchors are removed.


[[[BEGIN-SECTION:DESCRIPTION]]]
This is a command line tool provided with a Grackle installation. It fetches
and caches data files used by that associated Grackle version (the
``show-grackle-version`` subcommand shows that associated version).

In more detail:

- the tool caches downloaded data files in a central data-store that is
  accessible by pygrackle (and soon, Grackle, itself).

- within the data-store, this tool groups data-files inside of a directory
  named after the associated Grackle version. Consequently, the data-store
  (as a whole) is forwards AND backwards compatible, even if datafiles are
  renamed/mutated/added/deleted between Grackle versions.

- for the common scenario where the contents/names of data-files don't
  change between Grackle versions, this tool tries to deduplicate files.
  (This is accomplished with hard-links)

The location of the data is controlled by the ``GRACKLE_DATA_DIR``
environment variable. When this variable isn't specified, the tool uses
the operating-system recommendation for user-site-data. This location can
be queried with the ``show-data-dir`` subcomand.
[[[END-SECTION:DESCRIPTION]]]

[[[BEGIN-SECTION:MOTIVATION]]]
Why does this tool exist? Datafiles are required by **ANY** non-trivial
program (e.g. a simulation-code or python script) that invokes Grackle.

It is instructive to consider the historic experience of an end-user of one of
these programs. To build Grackle, they would typically clone the git repository
for Grackle (including the data files). To invoke their program, they would
manually specify the path to the downloaded data file. Frankly, this doesn't
seem so bad; the manual intervention is a minor inconvenience, at worst.
While it would be nice to eliminate the manual intervention, this it doesn't
seem it warrants development of a special tool.

Indeed, this is all true. Users who like this workflow can continue using it.
However, this manual management of datafiles becomes problematic in any
use-case that is marginally more complex. There are 3 considerations worth
highlighting:

1. **Portability:** Currently, there is no out-of-the-box approach for
   any program using Grackle configured to run on one computer to run on
   another machine without manual intervention.

   - If there are differences in how the machines are set up (e.g. where
     the data files are placed), the paths to the Grackle data file(s) need
     to be updated. This is relevant if you want to use a Pygrackle script
     on a different machine or if you want to use a configuration script to
     rerun a simulation (involving Grackle) on a different machine.

   - This is particularly noteworthy when it comes to automated testing! For
     example, before this tool existed, Pygrackle, made some assumptions that
     it was installed as an editable installation to run some examples. The
     test-suite of Enzo-E is another example where extra book-keeping is
     required for all test-problems that invoke Grackle.

2. **If the Grackle repository isn't present:** This includes the case where
   a user deletes the repository after installing Grackle. It is more
   important to consider the case where users are installing programs that use
   Grackle without downloading the repository (or, even if the repository is
   downloaded, it is done so without the user's knowledge). This latter case
   will become increasingly common as we make pygrackle easier to install.
   This is also plausible for cmake-builds of downstream projects that embed
   Grackle compilation as part of their build.

3. **Having multiple Grackle Versions Installed:** This is going to be
   increasingly common as Pygrackle becomes easier to install. Users have 2
   existing options in this case: (i) they maintain separate repositories of
   data files for each version or (ii) they assume that they can just use
   the newest version of the data-file repository. The latter option, has
   historically been true (and will probably continue to be true). But, it could
   conceivably lead to cases where people could unintentionally use a data-file
   created for a newer version of grackle. (While this likely won't be a
   problem, users should probably be explicitly aware that they are doing this
   on the off-chance that problems do arise).

This tool is a first step to addressing these cases.

Currently the tool just works for Pygrackle. There is an ongoing effort to add
functionality for the Grackle library, itself, to access the files managed by
this tool.
[[[END-SECTION:MOTIVATION]]]

[[[BEGIN-SECTION:INTERNALS-OVERVIEW]]]
We now turn our attention to describing how the organization of the data. It
all revolves around the idea of a **data-store** (more on that in a moment)

Protocol Version
++++++++++++++++
These internal details have an associated protocol-version, (queryable via the
``show-protocol-version`` subcommand). The details may change between protocol
versions. The protocol version will change very rarely (if it ever changes).

Data Directory
++++++++++++++
This is the directory that includes all grackle data. It's location is the
absolute path held by the ``GRACKLE_DATA_DIR`` environment variable. If the
variable isn't defined, it's placed at the location recommended by the
operating-system for user-site-data.

This contains several entries including the:

- a **user-data** directory. This directory currently isn't used yet, but
  it is reserved for users to put custom data-files in the future.

- the **data store** directory(ies). This is named
  ``data-store-v<PROTOCOL-VERSION>`` so that earlier versions of this
  tool will continue to function if we ever change the protocol.

Outside of the **user-data** directory, users should not modify/create
any files within Data Directory (it's okay to delete things).

Data Store
++++++++++
This is where we cache the data files managed by this system. This holds 1 or
more "version-directories". Each version-directory is named after a Grackle
version (**NOT** a Pygrackle version). A version-directory holds entries for
1 or more of the data-files that were shipped with the Grackle-version that
corresponds to the directory's name.

When this tool downloads new files, it applies basic deduplication logic:
After the file is downloaded, it checks whether a file already exists in an
adjacent version-directory with the same name and contents. If so, the
downloaded file is deleted and a hardlink is made to the pre-existing file.

- Users can be somewhat ignorant of our use of hardlinks. They are free to
  delete files (but the underlying data won't be removed until all of its links
  have been removed)

- This logic will **NOT** deduplicate renamed files. Doing that efficiently
  requires a different, more complex system (deleting the files is tricky). You
  can look at the commit history for an earlier version of this tool that can
  handle this scenario.
[[[END-SECTION:INTERNALS-OVERVIEW]]]

Sample Directory Structure
++++++++++++++++++++++++++
Down below, we sketch out what the directory-structure might look like:

[[[BEGIN:DIRECTORY-CARTOON]]]
${GRACKLE_DATA_DIR}/
 ├── data-store-v1/                      # <- the data-store
 │    ├── 3.3.1-dev/                     # <- a version-dir
 │    │    ├── CloudyData_UVB=FG2011.h5
 │    │    ├── ...
 │    │    └── cloudy_metals_2008_3D.h5
 │    └── 3.4.0/                         # <- another version-dir
 │         ├── CloudyData_UVB=FG2011.h5
 │         ├── ...
 │         └── cloudy_metals_2008_3D.h5
 └── user-data/                          # <- reserved for user data
      ├── ...
      └── ...
[[[END:DIRECTORY-CARTOON]]]

[[[BEGIN-SECTION:INTERNALS FILE REGISTRY]]]
The file registry refers to a small file that maps the names of Grackle's
datafiles with a checksum. At the time of writing we use sha256 checksums.

In the long term, the registry is important in 2 -ish contexts:

1. Embed the contents of the registry within the core Grackle library in order
   to provide the option to directly read in files from the data-store without
   being provided a full path. This embedding is handled by the build-system.
   PR #237 proposes an implmentation for this logic.

2. It's important for making make this tool work when:

   - this tool is part of the pygrackle module.
   - this tool is run as a portable single-file executable.
[[[END-SECTION:INTERNALS-REGISTRY]]]
"""

# For use as a standalone command line tool this MUST ONLY use python's built-in modules
import argparse
import contextlib
import filecmp
import functools
import hashlib
import io
import logging
from math import log10
import os
import re
import shutil
import stat
import sys
import traceback
from typing import IO, NamedTuple, Union
import urllib.request
from urllib.error import URLError, HTTPError


if sys.version_info < (3, 6, 1):  # 3.6.0 doesn't support all NamedTuple features
    raise RuntimeError("python 3.6.1 or newer is required")

_CHUNKSIZE = 8192  # default chunksize used for file operations

logger = logging.getLogger("grdata")
logger.setLevel(logging.DEBUG)
logger.addHandler(logging.StreamHandler())


class GrdataError(RuntimeError):
    pass


@contextlib.contextmanager
def _file_openner(f, mode, **kwargs):
    if isinstance(f, (str, bytes, os.PathLike)):
        with open(f, mode, **kwargs) as file:
            yield file
    else:
        yield f


def _progress_bar(tot_bytes, silent=False):
    """provides a function for drawing/updating progress bars"""
    ncols = shutil.get_terminal_size()[0] - 1
    power_div_3 = int(log10(tot_bytes) // 3) if tot_bytes > 0 else 0
    factor, unit = 1000.0**power_div_3, (" B", "KB", "MB", "GB")[power_div_3]
    # the output line has the form: '[<progress-bar>] <size>/<size> <unit>'
    fmt = "\r[{bar:{barlen}.{nfill}}] {size:.2f}" + f"/{tot_bytes / factor:.2f} {unit}"
    barlen = ncols - 19  # for context, 15 <= (len(fmt.format(...)) - barlen) <= 19
    suppress = (barlen < 1) or silent or not sys.stdout.isatty()
    bar = None if suppress else (barlen * "=")

    def _update(size):
        nonlocal bar
        if size is None and bar is not None:
            print(flush=True)
            bar = None
        elif bar is not None:
            nfill = int(barlen * (size / tot_bytes))
            val = fmt.format(bar=bar, barlen=barlen, nfill=nfill, size=size / factor)
            print(val, end="", flush=True)

    return _update


def _retrieve_url(url, dst, *, silent=False, chunksize=_CHUNKSIZE):
    """download the file from url to dst"""
    try:
        req = urllib.request.Request(url)
        with contextlib.ExitStack() as stack:
            out_file = stack.enter_context(open(dst, "wb"))
            response = stack.enter_context(urllib.request.urlopen(req))
            total_bytes = int(response.headers.get("Content-Length", -1))
            update_progress = _progress_bar(total_bytes, silent=silent)
            stack.callback(update_progress, size=None)

            # write downloaded data to a file
            downloaded_bytes = 0
            while True:
                update_progress(downloaded_bytes)
                block = response.read(chunksize)
                if not block:
                    break
                downloaded_bytes += len(block)
                out_file.write(block)
    except HTTPError as e:
        raise GrdataError(f"server can't fulfill request to fetch {url}: {e.code}")
    except URLError as e:
        raise GrdataError(f"server can't be reached to fetch {url}: {e.code}")


def _get_data_dir(system_str=None):
    """Returns a string specifying the data directory"""
    system_str = sys.platform if system_str is None else system_str
    if system_str.startswith("win32"):
        raise RuntimeError()

    # the following list specifies the order of directories in order of precedence
    # -> https://specifications.freedesktop.org/basedir-spec/latest/ instructs us to
    #    simply skip over XDG_DATA_HOME if it specifies a relative path
    # -> if you look back through the commit-history, you will see that an earlier
    #    version of the functionality used "$HOME/Library/Application Support/grackle"
    #    on macOS. For simplicity, we treat macOS like any other Unix
    data_dir_chain = [  # tuple fmt: (envvar, suffix_path, fatal_if_rel)
        ("GRACKLE_DATA_DIR", None, True),
        ("XDG_DATA_HOME", "grackle", False),
        ("HOME", ".local/share/grackle", True),
    ]

    for envvar, suffix_path, fatal_if_rel in data_dir_chain:
        # an earlier version of this function had special handling when envvar=="HOME"
        # & the variable wasn't defined (for simplicity we now treat that as an error)
        prefix = os.getenv(envvar, None)
        if prefix is None or len(prefix) == 0:
            continue
        elif prefix[0] == "~":
            msg = f"core grackle-lib forbids {envvar!r} path starting with '~'"
            raise RuntimeError(msg)
        elif prefix[0] != "/" and fatal_if_rel:
            raise RuntimeError(f"{envvar} holds a relative path")
        elif prefix[0] != "/":
            continue
        elif suffix_path is None:
            return os.path.normpath(prefix)
        else:
            return os.path.normpath(os.path.join(prefix, suffix_path))


class Config(NamedTuple):
    """Track basic configuration information."""

    data_dir: str
    base_url: str
    file_registry_file: Union[str, bytes, os.PathLike, IO]
    grackle_version: str
    protocol_version: str = "1"
    checksum_kind: str = "sha256"


def make_config_object(grackle_version, file_registry_file):
    """Construct the configuration object used for running the calculation

    Parameters
    ----------
    grackle_version : str
        the version of grackle (NOT pygrackle)
    file_registry_file : file or str or bytes or ``os.PathLike``
        Contains the file registry
    """
    _REPO_URL = "https://github.com/grackle-project/grackle_data_files"
    # this is the hash that holds the versions of the datafiles from the time when this
    # version of the file was shipped
    _CONTEMPORANEOUS_COMMIT_HASH = "9a63dbefeb1410483df0071eefcbff666f40816d"

    return Config(
        data_dir=_get_data_dir(),
        base_url=f"{_REPO_URL}/raw/{_CONTEMPORANEOUS_COMMIT_HASH}/input/",
        file_registry_file=file_registry_file,
        grackle_version=grackle_version,
    )


def _datastoredir_and_versiondir(conf):
    data_store_dir = os.path.join(conf.data_dir, f"data-store-v{conf.protocol_version}")
    version_dir = os.path.join(data_store_dir, conf.grackle_version)
    return data_store_dir, version_dir


def _setup_and_get_dirs(conf):
    """ensures file system is set up for fetching new files
    Returns dest-dir and list of directories to search for deduplication
    """

    # special handling for the fallback case on (non-macOS) Unix-like systems
    _path = os.path.normpath(os.path.expanduser(os.path.join("~", ".local", "share")))
    if os.fsdecode(conf.data_dir).startswith(_path):
        os.makedirs(_path, exist_ok=True)

    def _mkdir(path):
        if not os.path.isdir(path):
            os.mkdir(path)

    data_store_dir, dest_dir = _datastoredir_and_versiondir(conf)
    _mkdir(conf.data_dir)  # holds the data-store
    _mkdir(os.path.join(conf.data_dir, "user-data"))  # reserved for users' data
    _mkdir(data_store_dir)  # holds data for various Grackle versions
    _mkdir(dest_dir)  # holds data for current Grackle version

    with os.scandir(data_store_dir) as it:
        dirs = (entry.path for entry in it if entry.is_dir())
        search_dirs = [d for d in dirs if not os.path.samefile(dest_dir, d)]

    return dest_dir, search_dirs


def _parse_file_registry(f):
    """Read the file registry, as a dict from a text file

    Parameters
    ----------
    f : file or ``os.PathLike``
        Contains the data to be read in
    """
    # This format was choosen so that the contents could be directly injected into a
    # C file (and be interpreted as elements in an array of string-literals)
    with _file_openner(f, "r") as file:
        file_registry = {}
        for i, line in enumerate(file):  # iterater over lines
            line = line.rstrip()
            if (len(line) == 0) or line.isspace() or line.startswith("//"):
                continue
            elif line[0] != '"' or line[-2:] != '",':
                raise RuntimeError(f"Problem parsing {f}:{i + 1}\n   `{line}`")
            fname, cksum = line[1:-2].split()
            file_registry[fname] = cksum
    return file_registry


def calc_checksum(fname, alg_name, *, chunksize=_CHUNKSIZE):
    """Calculate the checksum for a given fname"""
    hash_calculator = hashlib.new(alg_name)
    with open(fname, "rb") as f:
        buffer = bytearray(chunksize)
        while True:
            nbytes = f.readinto(buffer)
            if nbytes == chunksize:
                hash_calculator.update(buffer)
            elif nbytes:  # equivalent to: (nbytes is not None) and (nbytes > 0)
                hash_calculator.update(buffer[:nbytes])
            else:
                break
    return ":".join([alg_name.lower(), hash_calculator.hexdigest()])


def matches_checksum(fname, checksum):
    alg_name, _ = checksum.split(":")
    return checksum == calc_checksum(fname, alg_name)


def _fetch(url, dst_dir, checksum, search_dirs):
    """Helper method to fetch a single file.

    Parameters
    ----------
    url: str
        The source url (e.g. file://path/to/fname, https://url/to/fname)
    dst_dir: str
        Path to the output directory
    checksum: str
        The checksum of the file.
    search_dirs: sequence of paths or None
        Should be ``None`` if ``dst_dir`` is an arbitrary location. Otherwise,
        this should point to a sequence of paths to all directories in the same
        datastore as dest_dir (other than ``dest_dir``). The other directories
        are used for deduplication. To disable deduplication, pass an empty
        sequence.
    """
    fname = os.path.basename(url)
    dst = os.path.join(dst_dir, fname)
    tmp_path = os.path.join(dst_dir, "_tempfile")
    if os.path.exists(tmp_path):
        os.remove(tmp_path)

    try:
        if url.startswith("file:"):
            shutil.copyfile(url[5:], dst=tmp_path)
        else:
            _retrieve_url(url, dst=tmp_path)
        if not matches_checksum(tmp_path, checksum):
            raise GrdataError(f"downloaded {fname} does't have expected checksum")
        os.rename(src=tmp_path, dst=dst)

        if search_dirs is not None:
            paths = list(
                filter(os.path.isfile, (os.path.join(d, fname) for d in search_dirs))
            )
            for path in paths:
                assert not os.path.samefile(path, dst)
                os.link(src=path, dst=tmp_path)
                if filecmp.cmp(tmp_path, dst, shallow=False):
                    os.remove(dst)
                    os.rename(src=tmp_path, dst=dst)
                    break
                os.remove(tmp_path)

            # set permissions to prevent accidental file-mutation/corruption (but still
            # allows deletion). We skip this for untracked data-dirs to avoid "are you
            # sure you want to delete this file?" prompts (on some platforms, like MacOS)
            os.chmod(dst, stat.S_IREAD | stat.S_IRGRP | stat.S_IROTH)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def fetch_all(config, *, fnames=None, override_dest=None):
    """Ensures that files in the specified registry are downloaded

    Parameters
    ----------
    config: Config
        basic configuration
    fnames: sequence, optional
        Optionally specifies a list of files, with corresponding
        registry entries to fetch. When this is ``None``, all files in
        the registry are fetched.
    override_dest: os.PathLike, optional
        When provided, fetched files are placed in this directory &
        no attempt is made to track the files as part of the standard
        grackle data-directory.
    """
    registry = _parse_file_registry(config.file_registry_file)

    # determine the files we will read (& their associated checksums)
    if fnames is None:
        fname_cksum_pairs = registry.items()
    else:
        for fname in filter(lambda fname: fname not in registry, fnames):
            raise ValueError(
                f"{fname} can't be downloaded since it doesn't have a registry entry."
                f"\n\nFiles with registry entries include: {list(registry.keys())!r}"
            )
        fname_cksum_pairs = ((fname, registry[fname]) for fname in fnames)

    # determine the output directory
    if override_dest is None:
        dest_dir, search_dirs = _setup_and_get_dirs(config)
        logger.info(f"destination dir: {dest_dir}")
    else:
        if not os.path.isdir(override_dest):
            os.mkdir(override_dest)
        dest_dir, search_dirs = override_dest, None
        logger.info("using override destination directory")

    logger.info(f"preparing to fetch files from: {config.base_url}")
    num_fetched = 0
    checksum_prefix = f"{config.checksum_kind}:"
    for fname, cksum_str in fname_cksum_pairs:
        dest_path = os.path.join(dest_dir, fname)

        if not cksum_str.startswith(checksum_prefix):
            if ":" in cksum_str:
                explanation = "was computed with a different algorithm"
            else:
                explanation = "is missing the `<cksum-alg>:` prefix"
            raise GrdataError(
                "To download a file, we must know the file's checksum computed with "
                f"the {config.checksum_kind} checksum algorithm. The provided checksum "
                f"for {fname} ({cksum_str}) {explanation}"
            )
        elif os.path.exists(dest_path):  # if the file already exists, we're done
            if not matches_checksum(dest_path, cksum_str):
                raise GrdataError(f"{dest_path} already exists but has wrong hash")
            continue
        else:
            logger.info(f"   fetching `{fname}`")
            url = os.path.join(config.base_url, fname)
            _fetch(url, dest_dir, cksum_str, search_dirs)
            num_fetched += 1

    if num_fetched == 0:
        logger.info("(no files needed to be retrieved)")


def _register_fetch_command(subparsers):
    parser_fetch = subparsers.add_parser(
        "fetch", help="fetch data files for associated Grackle version"
    )
    parser_fetch.add_argument(
        "fnames", nargs="*", help="Optional subset of files from the registry to fetch"
    )
    parser_fetch.add_argument(
        "--untracked-dest-dir",
        default=None,
        help=(
            "Download files to the specified arbitrary directory where we will not "
            "perform deduplication. This provided as a convenience. The specified "
            "directory should NOT be located inside the grackle data directory."
        ),
    )

    def fetch_command(args, config):
        fnames = None if len(args.fnames) == 0 else args.fnames
        fetch_all(config, fnames=fnames, override_dest=args.untracked_dest_dir)

    parser_fetch.set_defaults(func=fetch_command)


def calcreg_command(args, config):
    # calculate a new file registry for all files in the specifed directory
    try:
        with os.scandir(args.path) as scan_it:
            it = [(entry.name, entry.path) for entry in scan_it if entry.is_file()]
    except FileNotFoundError:
        raise ValueError(f"{args.path!r} doesn't specify a directory or file")
    except NotADirectoryError:
        it = [(os.path.basename(args.path), args.path)]

    lines = [f'"{name}  {calc_checksum(path, args.hash_alg)}",' for name, path in it]

    if not args.no_comments:
        print(f"""\
// This is a file registry generated by the grdata tool. To overwrite the file
// with an updated copy, install pygrackle and invoke:
//    python -m pygrackle calc-reg --hash-alg {args.hash_alg} {{dir}} > {{out-path}}
// where {{dir}} is the directory containing all files in the registry\
""")
    print(*lines, sep="\n")


def _register_calcreg_command(subparsers):
    descr = "compute & print file registry (file hash pairs) for a given directory."
    parser_calcreg = subparsers.add_parser("calc-reg", description=descr, help=descr)
    parser_calcreg.add_argument(
        "--no-comments", action="store_true", help="exclude comments from output"
    )
    _kinds = hashlib.algorithms_guaranteed
    parser_calcreg.add_argument(
        "--hash-alg", required=True, choices=_kinds, help="kind of checksum to compute"
    )
    parser_calcreg.add_argument(
        "path", help="path to the directory containing the files in the registry"
    )
    parser_calcreg.set_defaults(func=calcreg_command)


def _show_help(*args, **kwargs):
    # replace the anchors: they look like [[[<KIND>:<VALUE>]]]
    def repl(matchobj):
        anchor_kind, anchor_val = matchobj.group(1), matchobj.group(2)
        if anchor_kind == "BEGIN-SECTION":  # anchor_val is a section-name
            return "\n".join([anchor_val, len(anchor_val) * "-"])
        return "<SKIP>"

    _open, _close = re.escape("[[["), re.escape(r"]]]")
    p = re.compile(rf"^{_open}([-+0-9A-Za-z]+):([-,_+.! 0-9A-Za-z]+){_close}[ \t]*$")

    start = __doc__.index("[[[BEGIN-SECTION:DESCRIPTION]]]")
    it = (p.sub(repl, e).rstrip() for e in __doc__[start:].splitlines(keepends=False))

    for line in filter(lambda line: line != "<SKIP>", it):
        print(line)


def register_show_subcommand(name, short_descr, fn=None, *, subparsers):
    """register subcommand that simply shows information and then exits"""

    def _fn(args, config, attr):
        print(getattr(config, attr))

    func = functools.partial(_fn, attr=name.replace("-", "_")) if fn is None else fn
    _descr = f"show {short_descr} & exit"
    p = subparsers.add_parser(f"show-{name}", help=_descr, description=_descr)
    p.set_defaults(func=func)


def build_parser(prog_name):
    parser = argparse.ArgumentParser(
        prog=prog_name,
        description="Command line tool for downloading & caching Grackle's data files",
        epilog=f"Invoke `{prog_name} show-help` to get a detailed overview of the tool",
    )

    # Add flags for overriding some Config attrs (flags with the --testing-override-
    # prefix may be removed at any time in the future)
    def _add_override(attr, descr):
        prefix = "--testing-override-" if descr is None else "--override-"
        flag = prefix + attr.replace("_", "-")
        _help = argparse.SUPPRESS if descr is None else descr
        parser.add_argument(flag, dest=attr, default=argparse.SUPPRESS, help=_help)

    override_baseurl_descr = (
        "overrides the urlfiles are downloaded from. To copy files from a directory, "
        "you can specify `file:/path/to/dir`"
    )
    _add_override("base_url", override_baseurl_descr)
    _add_override("file_registry_file", None)
    _add_override("grackle_version", None)

    subparsers = parser.add_subparsers(required=True)
    _register_fetch_command(subparsers)
    _register_calcreg_command(subparsers)

    show = functools.partial(register_show_subcommand, subparsers=subparsers)
    show("grackle-version", "associated Grackle version")
    show("protocol-version", "associated data-store protocol version")
    show("checksum-kind", "associated checksum algorithm")
    show("data-dir", "path for data-dir (regardless of whether it exists)")
    show("help", "detailed help information", _show_help)

    def _showknownreg_command(args, config):
        with _file_openner(config.file_registry_file, "r") as f:
            print(*f.readlines(), sep="", end="")

    showknownreg_descr = "pre-registered file registry (for associated Grackle version)"
    show("known-reg", showknownreg_descr, _showknownreg_command)

    return parser


def main(default_config, prog_name, *, args=None):
    """Launch the command"""
    # modify the logger
    for handler in logger.handlers:
        logger.removeHandler(handler)
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(logging.Formatter("> %(message)s"))
    logger.addHandler(console_handler)

    # build parser & immediately parse arguments
    args = build_parser(prog_name).parse_args(args=args)

    # handle any overrides of config-attributes
    overrides = {f: getattr(args, f) for f in Config._fields if hasattr(args, f)}
    config = default_config._replace(**overrides)  # <- ok if len(overrides) == 0

    try:
        args.func(args, config=config)
        return 0
    except GrdataError as err:
        print(f"ERROR: {err.args[0]}")
        return 70  # https://www.man7.org/linux/man-pages/man3/sysexits.h.3head.html
    except BaseException:
        print("Unexpected error:", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        return 70  # https://www.man7.org/linux/man-pages/man3/sysexits.h.3head.html


# Here, we define machinery employed when used as a standalone program
# ====================================================================
def _exec_as_standalone_program():
    """executes the grdata program as a standlone program."""
    # To install the grdata as a standalone program, the build-system must:
    # - treat this file as a template & configure it with CMake's ``configure_file``
    #   command (or invoke ``configure_file.py`` under the classic build system) in
    #   order to substitute the names enclosed by the "at sign" symbol
    # - make resulting file executable
    # - install it into the bin directory alongside the grackle libraries

    variables = {  # holds the substituted variables
        "_GRDATA_GRACKLE_VERSION": "@_GRDATA_GRACKLE_VERSION@",
        "_GRDATA_FILE_REGISTRY_CONTENTS": """\
@_GRDATA_FILE_REGISTRY_CONTENTS@
""",
    }

    # check that substitution was sucessful
    for name, val in variables.items():
        # we use unicode escape sequence, \u0040, that python automatically converts
        # to the "at sign" to prevent the configure_file.py script (used by Grackle's
        # build-system) from falsely reporting an error
        if (name in val) or ("\u0040" in val) or (len(val) == 0) or (val.isspace()):
            raise RuntimeError(f"the build-system misconfigured the {name} variable")

    _CONFIG = make_config_object(
        grackle_version=variables["_GRDATA_GRACKLE_VERSION"],
        file_registry_file=io.StringIO(variables["_GRDATA_FILE_REGISTRY_CONTENTS"]),
    )
    sys.exit(main(default_config=_CONFIG, prog_name="grdata"))


if __name__ == "__main__":
    _exec_as_standalone_program()
