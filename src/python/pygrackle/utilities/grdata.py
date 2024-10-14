#!/usr/bin/env python3

# A tool for managing grackle data files. (More details provided after imports)
#
# This file should be usable as both (i) a part of pygrackle and (ii) a standalone
# command line tool (when pygrackle IS NOT installed)
#
# To support scenario 1, this CAN ONLY use python's built in modules.

import argparse
import contextlib
import filecmp
import hashlib
import io
from math import log10
import os
import re
import shutil
import stat
import sys
import traceback
from typing import IO, NamedTuple, Optional, Union
import urllib.request
from urllib.error import URLError, HTTPError
import warnings

# Down below, we provide a detailed description that serves 3 purposes
#   1. to act as a description of this files contents for developers
#   2. to serve as documentation on the website
#   3. to serve as queryable documentation via the `help` subcommand
#
# The text enclosed in triple braces serves 2 purposes:
# -> it is designed to be anchors used by sphinx to include the documentation.
# -> while executing the `help` subcommand, anchors of the format
#    `[[[BEGIN-SECTION:<name>]]]` will be replaced with a section title `<name>`,
#    and all other anchors are removed.

_EXTENDED_DESCRIPTION = """\
[[[BEGIN-SECTION:DESCRIPTION]]]
This is a management system for managing Grackle data files. The command line
interface provides commands to fetch these data files, list all of the
available data, and delete the data.

The system stores the data files at a single global location. (Grackle,
itself, can access files from this location).

The key feature of this system is its support for versioning:

- it is able to support management of sets of datafiles (associated with
  different grackle versions) where the datafiles have been renamed,
  modified, or deleted between Grackle versions.

- additionally, the system implements deduplication for the (very common)
  scenario when the contents of a file are unchanged between grackle
  versions.

One minor caveat: a given version of this tool is ONLY able to download
data for the grackle version specified by the ``--version-grackle`` flag
(i.e. this is the grackle version that the tool ships with). However, it
does support listing and deleting data associated with other grackle
versions.

The location of the data is controlled by the ``GRACKLE_DATA_DIR``
environment variable. When this variable isn't specified, the tool uses
the operating-system recommendation for user-site-data. This location can
be queried with the ``getpath`` subcomand.
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
We now turn our attention to describing how the internals of the
management system work.

Fundamentally, the data management system manages a **data store**.
We will return to that in a moment.

Protocol Version
++++++++++++++++

This internal logic has an associated protocol-version, (you can query
this via the ``--version-protocol`` flag). The logic may change between
protocol versions. The protocol version will change very rarely (if it
ever changes at all)

Data Directory
++++++++++++++

This is simply the data directory that includes all grackle data. This path
is given by the ``GRACKLE_DATA_DIR`` environment variable, if it exists.
Otherwise it defaults to the operating-system's recommendation for
user-site-data.

This contains several entries including the:

  - a **user-data** directory. This directory currently isn't used yet, but
    it is reserved for users to put custom data-files in the future.

  - a **tmp** directory (used by the data-management tool)

  - it sometimes holds a lockfile (used to ensure that multiple instances of
    this tool aren't running at once)

  - the **data store** directory(ies). This is named
    ``data-store-v<PROTOCOL-VERSION>`` so that earlier versions of this
    tool will continue to function if we ever change the protocol. (Each of
    these directories are completely independent of each other).

Outside of the **user-data** directory, users should not modify/create/delete
any files within Data Directory (unless the tool instructs them to).

Data Store
++++++++++

This is where we track the data files managed by this system. This holds a
directory called **object-store** and 1 or more "version-directories".

The primary-representation of each file is tracked within the ``object-store``
subdirectory.

- The name of each item in this directory is a unique key. This key is the
  file’s SHA-1 checksum.

- Git internally tracks objects in a very similar way (they have historically
  used SHA-1 checksums as unique keys). The chance of an accidental collision
  in the checksum in a large Git repository is extremely tiny. It was only 10 or
  12 years after Git was created that the developers started worrying about
  collisions (and they are primarily concerned with intentional collisions from
  maclicious actors).

Each version-directory is named after a Grackle version (**NOT** a Pygrackle
version).

- a given version directory holds data-file references.
- the references have the contemporaneous name of each of the data-files that
  was shipped with the Grackle-version that corresponds to the directory's name.
- each reference is linked to the corresponding file in the ``object-store``.

When a program outside of this tool accesses a data-file, they will **ONLY**
access the references in the version-directory that shares its name with the
version of Grackle that the program is linked against.

This tool makes use of references and the ``object-store`` to effectively
deduplicate data. Whenever this tool deletes a "data-file" reference it will
also delete the corresponding file from the ``object-store`` if it had no other
references. We choose to implement references as "hard links" in order to make
it easy to determine when a file in ``object-store`` has no reference.
[[[END-SECTION:INTERNALS-OVERVIEW]]]

Sample Directory Structure
++++++++++++++++++++++++++

Down below, we sketch out what the directory-structure might look like:

[[[BEGIN:DIRECTORY-CARTOON]]]
GRACKLE_DATA_DIR/
 ├── data-store-v1/                      # <- the data-store
 │    ├── 3.3.1-dev/                     # <- a version-dir
 │    │    ├── CloudyData_UVB=FG2011.h5
 │    │    ├── ...
 │    │    └── cloudy_metals_2008_3D.h5
 │    ├── 3.4.0/                         # <- another version-dir
 │    │    ├── CloudyData_UVB=FG2011.h5
 │    │    ├── ...
 │    │    └── cloudy_metals_2008_3D.h5
 │    └── object-store/                  # <- the object-store
 │         ├── ...
 │         └── ...
 ├── tmp/                                # <- reserved for scratch-space
 ├── user-data/                          # <- reserved for user data
 │    ├── ...
 │    └── ...
 └── lockfile                            # <- temporary file
[[[END:DIRECTORY-CARTOON]]]
"""


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


# define some constants
# =====================

# default chunksize used for file operations
_CHUNKSIZE = 4096

# the name of the subdirectory in a data-store where we handle deduplication
_OBJECT_STORE_SUBDIR = "object-store"

# Alternative to `None` for specifying that a value wasn't specified. This is primarily
# used as a default value for an optional command line flag that requires an argument.
# In that case, a value of _UNSPECIFIED means that the flag wasn't specified while a
# value of `None` means that the flag doesn't have an associated value.
_UNSPECIFIED = object()


# Version check and define some backports
# =======================================

if sys.version_info < (3, 6, 1):
    # 3.6.0 doesn't support all NamedTuple features
    raise RuntimeError("python 3.6.1 or newer is required")

elif sys.version_info < (3, 7, 0):

    class nullcontext:
        def __init__(self, enter_result=None):
            self.enter_result = enter_result

        def __enter__(self):
            return self.enter_result

        def __exit__(self, *args):
            pass

else:
    nullcontext = contextlib.nullcontext


class GenericToolError(RuntimeError):
    pass


class ToolConfig(NamedTuple):
    """Tracks basic information about this tool"""

    grackle_version: str
    protocol_version: str = "1"
    checksum_kind: str = "sha1"


def _ensure_all_removed(fnames):
    for fname in fnames:
        try:
            os.remove(fname)
        except FileNotFoundError:
            continue


_MAX_BAR = 160 * "="


@contextlib.contextmanager
def _progress_bar(ncols, total_bytes, *, use_dummy=False):
    """
    ContextManager that provides a function used for drawing/updating progress bars

    If the program wasn't excuted from a shell, or the caller want to draw
    too few columns, the returned function does nothing.
    """
    # the main template is '[<progress-bar>] <size>/<size> <unit>'
    # -> <progress-bar> is some fraction of _FULL_BAR and empty space
    # -> <size> describes the current/total download size (takes up to 6 characters)
    # -> <unit> is 1 or 2 characters
    # -> thus, we need 19 characters for everything other than <progress-bar>
    bar_len = min(len(_MAX_BAR), ncols - 19)

    if use_dummy or (total_bytes <= 0) or (bar_len <= 0) or not sys.stdout.isatty():
        use_dummy = True

        def _update(size):
            return None
    else:
        power_div_3 = int(log10(total_bytes)) // 3
        factor, unit = 1000.0**power_div_3, ("B", "KB", "MB", "GB")[power_div_3]
        fmt = "\r[{bar:{len}.{nfill}}] {size:.2f}" + f"/{total_bytes/factor:.2f} {unit}"

        def _update(size):
            nfill = bar_len * int(size / total_bytes)
            val = fmt.format(bar=_MAX_BAR, len=bar_len, nfill=nfill, size=size / factor)
            print(val, end="", flush=True)

    try:
        yield _update
    finally:
        # always execute this clause when exiting the context manager. If an exception
        # caused the exit, it will be re-raised after this clause
        if not use_dummy:
            print(flush=True)


def _retrieve_url(url, dest, fname, *, use_progress_bar=True, chunksize=_CHUNKSIZE):
    """
    download the file from url to dest

    Note
    ----
    Online discussion about calling `response.read(chunksize)`, where
    `response` is the context manager object produced by
    `url.request.urlopen`, seems to strongly imply that this limits the
    amount of data read from the http request into memory at a given
    point in time. However, the documentation seems vague on this point.

    This is unlikely to ever be a problem (the biggest file we need is
    currently under 10 Megabytes). However, if it does become a problem,
    we have 2 options:
      1. we could conditionally fall back to ``curl`` (cli tool) or
         ``Requests`` (python package) if they are present on the system
      2. we could craft custom http requests
    """
    ncols = shutil.get_terminal_size()[0] - 1
    req = urllib.request.Request(url)
    try:
        with contextlib.ExitStack() as stack:
            # enter context managers for http-response, progress-bar, & output-file
            response = stack.enter_context(urllib.request.urlopen(req))
            total_bytes = int(response.headers.get("Content-Length", -1))
            update_progress = stack.enter_context(
                _progress_bar(ncols, total_bytes, use_dummy=not use_progress_bar)
            )
            out_file = stack.enter_context(open(dest, "wb"))

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
        raise GenericToolError(
            f"The server couldn't fulfill the request for retrieving {fname} from "
            f"{url}.\nError code: {e.code}"
        )
    except URLError as e:
        raise GenericToolError(
            f"The server couldn't be reached while trying to retrieve {fname} from "
            f"{url}.\nError code: {e.code}"
        )


class Fetcher(NamedTuple):
    """Encodes information for fetching data files

    Note
    ----
    Right now, we always assume that we want to support downloading from
    GitHub, but in the future, we can also support fetching from a
    directory
    """

    base_path: str
    holds_url: bool

    @classmethod
    def configure_GitHub_url(cls, data_repository_url, contemporaneous_git_hash):
        repo_url = data_repository_url
        repo_version = contemporaneous_git_hash
        # we could also use the name of a branch (instead of a commit-hash) if we think
        # that would be better
        return cls(base_path=f"{repo_url}/raw/{repo_version}/input/", holds_url=True)

    @classmethod
    def configure_src_dir(cls, dir_path):
        return cls(base_path=dir_path, holds_url=False)

    def __call__(self, fname, checksum, checksum_kind, dest_dir):
        """
        Retrieve the file named ``fname`` to a location dest_dir

        Returns
        -------
        full_path: str
            Upon success, we return the full path of the newly fetched file

        Notes
        -----
        We follow the following procedure (inspired by pooch):
          1. downloads the file to a temporary location
          2. verifies that the checksum is correct
          3. move the file to the appropriate destination

        This provides robust behavior if the program is interupted. In
        principle, we could combine steps 1 and 2. But, there may be some
        minor benefits to keeping our procedure like this (theoretically,
        we may be more likely to catch corruption from harddrive hardware
        errors).
        """
        src = os.path.join(self.base_path, fname)
        dst = os.path.join(dest_dir, fname)
        _pretty_log(f"-> fetching `{fname}`", indent_all=True)
        tmp_name = os.path.join(dest_dir, "_tempfile")
        # tmp_name can safely be removed if it exists (it only exists if this logic
        # previously crashed or was interupted by SIGKILL)
        _ensure_all_removed([tmp_name])

        try:
            if self.holds_url:
                _retrieve_url(src, tmp_name, fname)
            else:
                # copy the file
                shutil.copyfile(src, tmp_name)
            if not matches_checksum(tmp_name, checksum_kind, checksum):
                if matches_checksum(src, checksum_kind, checksum):
                    raise GenericToolError(
                        f"while copying from {src}, data may have been corrupted"
                    )
                raise GenericToolError(f"{src} does't have the expected checksum")
            os.rename(tmp_name, dst)

        finally:
            _ensure_all_removed([tmp_name])


class DataStoreConfig(NamedTuple):
    """Track basic configuration information

    In principle, this information is intended to be a little more
    flexible and might not be known as early as ToolConfig.
    """

    data_dir: str
    store_location: str
    checksum_kind: str
    default_fetcher: Fetcher
    file_registry_file: Union[str, bytes, os.PathLike, IO, None]

    @property
    def tmp_dir(self):
        """Used for hardlink test and scratch-space"""
        return os.path.join(self.data_dir, "tmp")

    @property
    def user_data_dir(self):
        """Reserved for user data"""
        return os.path.join(self.data_dir, "user-data")


def _get_platform_data_dir(appname="grackle", system_str=None):
    """Returns a string specifying the default data directory

    All of these choices are inspired by the API description of the platformdirs python
    package
        * we only looked at online documentation:
          https://platformdirs.readthedocs.io/en/latest/
        * we have NOT read any source code
    """
    if system_str is None:
        system_str = sys.platform
    if system_str.startswith("win32"):
        raise RuntimeError()
    elif system_str.startswith("darwin"):
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
        if env_str[0] == "~":
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
        if not manual_choice[:2] != "~/":  # for parity with C-version of this function
            raise RuntimeError(
                "can't expand can't expand env-variable, GRACKLE_DATA_DIR when "
                "it starts with `~user/` or just contains `~`"
            )
        return os.path.expanduser(manual_choice)
    else:
        return manual_choice


@contextlib.contextmanager
def _file_openner(f, mode, **kwargs):
    """Open a file or pass through an already open file"""
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
    if path is None:
        yield f
    else:
        with open(path, mode, **kwargs) as fobj:
            yield fobj


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

    with _file_openner(f, "r") as file:
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


class LockFileExistsError(FileExistsError):
    pass


class LockFileContext:
    """Reentrant context manager that creates a "lockfile".

    The context-manager will delete the file when we finish. If the lock
    already exists, the program will abort with an explanatory error
    (this ensures that only 1 copy of the program will try to run at a
    time).

    Examples
    --------
    To use this you might invoke:

    >>> dir_lock = LockFileContext("path/to/lockfile")
    >>> with dir_lock:
    ...     # do something critical

    This is reentrant in the sense that you can perform something like the
    following (the real value here is that you can mover internal
    with-statement inside of functions)

    >>> dir_lock = LockFileContext("path/to/lockfile")
    >>> with dir_lock:
    ...     # do something critical
    ...     with dir_lock:
    ...         # do something else critical
    """

    def __init__(self, lock_file_path):
        self.lock_file_path = lock_file_path

        # the following is always non-negative. It can exceed 1 if the same context
        # manager is used in nested with-statements
        self._acquisition_count = 0

    def locked(self):
        return self._acquisition_count > 0

    def __enter__(self):
        if self._acquisition_count == 0:
            # try to acquire the lock (by trying to create the file)
            try:
                f = open(self.lock_file_path, "x")
                f.close()
            except FileExistsError as err:
                raise LockFileExistsError(
                    err.errno,
                    err.strerror,
                    err.filename,
                    getattr(err, "winerror", None),
                    err.filename2,
                ) from None
        else:
            # this is a nested with-statement, in a process that already owns the lock
            pass

        self._acquisition_count += 1  # only executed if process owns the lock

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type is FileExistsError:
            return False
        elif self._acquisition_count <= 0:
            raise RuntimeError("the contextmanager has a totally invalid state!")
        elif self._acquisition_count == 1:
            os.remove(self.lock_file_path)
        self._acquisition_count -= 1
        return False  # if an exception triggered the exitting of a context manager,
        # don't suppress it!


def standard_lockfile(data_store_config):
    return LockFileContext(os.path.join(data_store_config.data_dir, "lockfile"))


def calc_checksum(fname, alg_name, *, chunksize=_CHUNKSIZE):
    """Calculate the checksum for a given fname"""
    # construct the object to track intermediate state of the checksum
    # calculation as we stream through the data
    hash_obj = hashlib.new(alg_name)
    with _file_openner(fname, "rb") as f:
        if f is fname:
            f.seek(0, os.SEEK_SET)

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


def _pretty_log(arg, *, indent_all=False):
    """indent messages so it's clear when multiline messages are a single thought"""
    lines = arg.splitlines()
    if len(lines) and not indent_all:
        formatted = [f"-- {lines[0]}"] + [f"   {e}" for e in lines[1:]]
    else:
        formatted = [f"   {e}" for e in lines]
    print(*formatted, sep="\n")


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
    def is_supported(dirname):
        """returns whether the OS (and filesystem supports hardlinks)"""

        fnames = [os.path.join(dirname, f"linktest_f{i}.txt") for i in [0, 1]]
        _ensure_all_removed(fnames)

        try:
            _contents = "THIS IS SOME TEST DATA"
            with open(fnames[0], "w") as f:
                f.write(_contents)
            os.link(fnames[0], fnames[1])
            os.remove(fnames[0])
            with open(fnames[1], "r") as f:
                support_hardlinks = f.read() == _contents
        except OSError:
            support_hardlinks = False
        finally:
            _ensure_all_removed(fnames)
        return support_hardlinks

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


def _ensure_data_dir_exists(data_store_config):
    """Creates the data_dir if it doesn't exist

    the data_dir is a directory that contains:
     -> the data-store directory for data managed by the current protocol version
     -> (possibly) data-store directories for data managed by other protocol version
     -> (possibly) a directory called `user-data/` where users can put custom data
    """
    _ensure_exists(data_store_config.data_dir, "that will hold all Grackle data")

    # even though it isn't used for anything right now, make the directory that is
    # reserved for user content
    _ensure_exists(data_store_config.user_data_dir, "reserved for user-defined data")

    # primarily for testing whether hard-links are supported
    _ensure_exists(data_store_config.tmp_dir, "reserved for scratch-space")


def get_version_dir(tool_config, data_store_config):
    return os.path.join(data_store_config.store_location, tool_config.grackle_version)


def get_object_dir(data_store_config):
    return os.path.join(data_store_config.store_location, _OBJECT_STORE_SUBDIR)


class VersionDataManager(NamedTuple):
    """
    Actually manages downloads of files to a directory where the
    directory and files are associated with single Grackle version.

    Instances of this class support 2 modes of operation:
      1. The instance manages a data file directory that is part of the
         larger data-management system. This data-management system may
         include multiple version directories and data-files are
         deduplicated.
      2. The instance manages a data file directory that is completely
         isolated from any data-management system (i.e. there is no
         deduplication)

    The first mode is the primary usecase of this class. (The second mode
    is mostly provided as a convenience)

    Warnings
    --------
    This should not be considered part of a public API. The names and
    existence of all attributes and methods are subject to change

    Notes
    -----
    A major motivating factor in the design was providing the capacity
    to create the necessary directories only when absolutely necessary
    (i.e. when we are about to download data)

    Some future methods that might be worth implmenting

      * a method to download a single file

      * a method to check the validity of a single Version file (i.e. it
        ONLY contains files listed in the specified registry, all files
        match the specified checksum, AND they are all properly linked to
        a file in the object directory)
    """

    # Path to output directory, where the file-name matches the name given in the
    # registry and is known by the associated grackle-version.
    version_dir: str

    # data_store_config holds a little more information than we actually need
    # -> we may chooise to redefine this in the future
    data_store_config: Optional[DataStoreConfig]

    # encodes the configuration (and logic) for fetching the files
    fetcher: Fetcher

    @classmethod
    def create(
        cls,
        tool_config,
        data_store_config,
        *,
        untracked_dest_dir=None,
        override_fetcher=None,
    ):
        """
        Create a new instance

        Parameters
        ----------
        tool_config : ToolConfig
        data_store_config : DataStoreConfig
        untracked_dest_dir : str, optional
            When specified the constructed object will be used to manage
            fetched files are placed in the specified
            directory and no attempt is made to track the file as part of
            the data-directory.
        override_fetcher : Fetcher, optional
            When specified, this fetcher is used in place of the standard
            default fetcher provided by data_store_config.
        """

        fetcher = override_fetcher
        if fetcher is None:
            fetcher = data_store_config.default_fetcher

        if untracked_dest_dir is None:
            version_dir = get_version_dir(tool_config, data_store_config)
        else:
            version_dir = untracked_dest_dir
            data_store_config = None

        return cls(
            version_dir=version_dir,
            data_store_config=data_store_config,
            fetcher=fetcher,
        )

    def manages_untracked_data(self):
        return self.data_store_config is None

    def _object_dir(self):
        """
        If the instance manages data as part of a larger data management
        system, then this method returns the path to the object directory.
        This is the directory where checksum names are used as filenames.
        (Linking with files in this directory is the mechanism used to aid
        deduplication)
        """
        if self.data_store_config is None:
            return None
        return get_object_dir(self.data_store_config)

    def _setup_file_system(self):
        """
        helper function that ensures that the file system is set up for
        fetching new files and returns the configured lockfile context
        manager (it isn't locked yet)
        """
        if self.manages_untracked_data():
            lockfile_ctx = nullcontext()
        else:
            _ensure_data_dir_exists(self.data_store_config)

            lockfile_ctx = standard_lockfile(self.data_store_config)
            with lockfile_ctx:
                # let's validate we can actually use hardlinks
                if not hasattr(os, "link"):
                    raise GenericToolError("operating system doesn't support hardlinks")
                elif not _HardlinkStrat.is_supported(self.data_store_config.tmp_dir):
                    raise GenericToolError("The file system does not support hardlinks")

                # a little more set up
                _ensure_exists(
                    self.data_store_config.store_location, "that holds the data-store"
                )
                _ensure_exists(self._object_dir(), "")

            assert not lockfile_ctx.locked()  # sanity check!

        with lockfile_ctx:
            _ensure_exists(
                self.version_dir, "that holds data for current Grackle version"
            )

        return lockfile_ctx

    def _fetch_file(self, fname, full_checksum_str, *, lockfile_ctx=None):
        """
        Helper method to fetch a single file. Provides the full path

        Parameters
        ----------
        fname : str
            The name of the file to be fetched
        full_checksum_str : str
            The checksum of the file.
        lockfile_ctx : LockFileContext, optional
            When this is None, the calling process doesn't already own the
            lock for the data-directory (and this function will try to
            acquire the lock)

        Returns
        -------
        any_work : bool
            Indicates whether any work was actually required to fetch the
            file. the file was freshly fetched. ``True`` indicates that we
            actually needed to go get the file, while ``False`` denotes
            that the file already existed
        full_fname : str
            Specifies the absolute path to the file. When a tracked file
            is fetched, then this is path to the file entry where
            `os.path.basename(full_fname)` is equal `fname`
        """

        if lockfile_ctx is None:
            lockfile_ctx = self._setup_file_system()

        # extract the checksum_kind and string that are stored in the registry
        if ":" in full_checksum_str:
            cur_cksum_kind, checksum = full_checksum_str.split(":")
        else:
            raise ValueError(
                f"the checksum for {fname} does not specify the checksum kind"
            )

        # when tracking files as part of the data file management system, there are
        # strict requirements on the kind of checksum that is used.
        # -> This is because the entries in the object directory (used for
        #    deduplication) use the checksum string (without the algorithm tag) as
        #    file names.
        # -> here, we check this requirement
        if self.manages_untracked_data():
            req_cksum_kind = None
        else:
            req_cksum_kind = self.data_store_config.checksum_kind
            if cur_cksum_kind != req_cksum_kind:
                raise ValueError(
                    "To download a file as part of Grackle's data file management "
                    "system, we must know the file's checksum computed with the "
                    f"{req_cksum_kind} checksum algorithm. The provided checksum "
                    f"for {fname} was computed with the {cur_cksum_kind} algorithm."
                )

        # now we actually fetch the file (if necessary)
        with lockfile_ctx:
            # get the full path to the downloaded file
            full_fname = os.path.join(self.version_dir, fname)

            # if the file already exists we are done
            if os.path.exists(full_fname):
                if not matches_checksum(full_fname, cur_cksum_kind, checksum):
                    raise RuntimeError(
                        f"{full_fname} already exists but has the wrong hash"
                    )
                # when we're handling tracked data, we could theoretically check whether
                # this data is properly linked. But I don't think I would want to take
                # any kind of action if it isn't properly tracked (the most action I
                # would recommend is providing the user with a warning)
                return (False, full_fname)

            # download the file
            fetcher = self.fetcher
            fetcher(
                fname,
                checksum=checksum,
                checksum_kind=cur_cksum_kind,
                dest_dir=self.version_dir,
            )
            if not self.manages_untracked_data():
                # by changing permissions, certain platforms (like MacOS) will ask
                # users "are you sure you want to delete this file" when the ``rm``
                # command is used:
                # -> this is desirable for files managed by the grackle data management
                #    systems (when we don't want users removing individual files)
                # -> but, we avoid doing this for untracked data files (because it's
                #    just an annoyance for the end-user)
                os.chmod(full_fname, _IMMUTABLE_MODE)

            if not self.manages_untracked_data():
                # handle deduplication (since the data files are part of the larger
                # grackle data file management system)
                cksum_fname = os.path.join(self._object_dir(), checksum)

                try:
                    _HardlinkStrat.deduplicate(full_fname, cksum_fname)

                    # not strictly necessary, but doing this for safety reasons
                    os.chmod(cksum_fname, _IMMUTABLE_MODE)

                except Exception as err:
                    # remove full_fname
                    # -> we don't want users to use it before resolving the issues
                    # -> We also want to make the errors as reproducible as possible
                    #    (ideally, rerunning the command should produce the same error
                    #    if you haven't changed anything)
                    os.remove(full_fname)
                    if os.path.is_file(cksum_fname) and not isinstance(err, ValueError):
                        raise err

                    # this should only happens when full_fname and cksum_fname both
                    # exist, but aren't perfect matches of each other. Here, we try to
                    # provide a more informative error message
                    if not matches_checksum(cksum_fname, req_cksum_kind, checksum):
                        raise GenericToolError(f"""\
A file (used for deduplication) that already existed on disk
   `{cksum_fname}`
which is probably a version of `{fname}`,
doesn't have the appropriate {req_cksum_kind} checksum.
-> expected: {calc_checksum(cksum_fname, req_cksum_kind)}
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
        return (True, full_fname)

    def fetch_all(self, registry, *, fnames=None):
        """
        Ensures that files in the specified registry are downloaded

        Parameters
        ----------
        registry : dict
            maps file names to associated checksums
        fnames : sequence, optional
            Optionally specifies a list of files, with corresponding
            registry entries to fetch. When this is ``None``, all files in
            the registry are fetched.
        """

        if fnames is None:
            fname_cksum_pairs = registry.items()
        else:
            for fname in filter(lambda fname: fname not in registry, fnames):
                raise ValueError(
                    f"{fname} is not the name of a file with a registry "
                    "entry. Thus it can't be downloaded.\n\nFiles with "
                    f"registry entries include: {list(registry.keys())!r}"
                )
            fname_cksum_pairs = ((fname, registry[fname]) for fname in fnames)

        # ensure all needed directories exist and fetch the lockfile context manager
        lockfile_ctx = self._setup_file_system()

        with lockfile_ctx:
            num_fetched = 0
            _pretty_log(f"preparing to fetch files from: {self.fetcher.base_path}")
            for fname, full_checksum_str in fname_cksum_pairs:
                any_work, _ = self._fetch_file(
                    fname, full_checksum_str, lockfile_ctx=lockfile_ctx
                )
                num_fetched += any_work

        if num_fetched == 0:
            _pretty_log("-> no files needed to be retrieved", indent_all=True)


def fetch_command(args, tool_config, data_store_config):
    override_fetcher = None
    if args.from_dir is not None:
        override_fetcher = Fetcher.configure_src_dir(args.from_dir)

    fnames = None if len(args.fnames) == 0 else args.fnames

    man = VersionDataManager.create(
        tool_config=tool_config,
        data_store_config=data_store_config,
        untracked_dest_dir=args.untracked_dest_dir,
        override_fetcher=override_fetcher,
    )
    registry = _parse_file_registry(data_store_config.file_registry_file)
    man.fetch_all(registry, fnames=fnames)


def _register_fetch_command(subparsers):
    parser_fetch = subparsers.add_parser(
        "fetch",
        help=(
            "fetch data files if we don't already have the data for the "
            "associated version of grackle"
        ),
    )
    parser_fetch.add_argument(
        "fnames",
        nargs="*",
        help=(
            "Optionally specify the names of files that should be fetched. "
            "Each listed file must have a corresponding entry in the file "
            "registry used by this tool. If no files are specified, then the "
            "tool will fetch every known file."
        ),
    )
    parser_fetch.add_argument(
        "--untracked-dest-dir",
        default=None,
        help=(
            "This flag can be used to instruct to download files to an "
            "arbitrary directory where the data files will NOT be stored as "
            "part of the data file management system. This provided as a "
            "convenience. The specified directory should NOT be located "
            "inside the grackle data directory. The tool also won't use any "
            "kind of lock files (thus, it's the user's responsibility to "
            "ensure that only a single process is modifying the specified "
            "directory at any given time. This is provided mostly as a "
            "convenience. Management of the files in this directory is "
            "outside this tool's scope."
        ),
    )
    parser_fetch.add_argument(
        "--from-dir",
        default=None,
        help=(
            "optionally specify a path to a directory where we copy the files from "
            "(instead of downloading them)"
        ),
    )
    parser_fetch.set_defaults(func=fetch_command)


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

    def always_true(*args):
        return True

    if ftype is None:
        has_ftype = always_true
    elif ftype == "dir":
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


def rm_command(args, tool_config, data_store_config):
    """Logic for removing files"""
    if args.vdata is _UNSPECIFIED:
        # this means that we are removing the whole data store
        if not args.data_store:
            raise RuntimeError("SOMETHING WENT HORRIBLY, HORRIBLY WRONG")

        _descr = os.path.basename(data_store_config.store_location)
        target_path = data_store_config.store_location
        operation_description = (
            f"deleting ALL files in the data-store associated with this tool, {_descr}"
        )
        if not os.path.isdir(target_path):
            raise GenericToolError(
                "intended to recursively delete all contents of the associated "
                "data-store. But no such directory can be found."
            )

        fn = shutil.rmtree

    else:
        if args.vdata is None:
            target = tool_config.grackle_version
            _descr = f"associated with this tool (`{tool_config.grackle_version}`)"
        else:
            target = args.vdata
            _descr = f"`{target}`"
        target_path = os.path.join(data_store_config.store_location, target)
        operation_description = (
            f"deleting all data file references for the grackle-version {_descr}. "
            "Any files for which the reference-count drops to zero will also be "
            "removed."
        )

        if not os.path.isdir(target_path):
            raise GenericToolError(
                "intended to delete all data-file references for the grackle-version "
                f"{_descr}, but no such data is tracked in the data-store."
            )

        def fn(path):
            object_dir = os.path.join(
                data_store_config.store_location, _OBJECT_STORE_SUBDIR
            )
            if not os.path.isdir(object_dir):
                raise RuntimeError(
                    "SOMETHING IS HORRIBLY WRONG!!! THE {object_dir} IS MISSING"
                )

            # we throw an err if this directory contains some unexpected stuff
            it = direntry_iter(path, ftype="file", mismatch="eager_err")
            for name, full_path in it:
                # get path to corresponding hardlinked file in _OBJECT_STORE_SUBDIR
                checksum = calc_checksum(full_path, alg_name=tool_config.checksum_kind)
                cksum_fname = os.path.join(object_dir, checksum)
                cksum_fname_exists = os.path.isfile(cksum_fname)

                if not cksum_fname_exists:
                    warnings.warn(
                        "Something weird has happened. There is no deduplication file "
                        f"associated with {full_path}"
                    )
                os.remove(full_path)
                if cksum_fname_exists:
                    _HardlinkStrat.remove_if_norefs(cksum_fname)
            os.rmdir(path)

    with standard_lockfile(data_store_config):
        if not args.force:
            _pretty_log(
                f"{operation_description}\n"
                "-> essentially, we are recursively removing\n"
                f"     `{target_path}`\n"
                "-> to actually perform this command, pass the --force flag"
            )
        else:
            fn(target_path)


def _register_rm_command(subparsers):
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


def lsversions_command(args, tool_config, data_store_config):
    if not os.path.exists(data_store_config.store_location):
        print("there is no data")
    with standard_lockfile(data_store_config):
        it = direntry_iter(
            data_store_config.store_location,
            ftype="dir",
            mismatch="lazy_err",
            ignore=[_OBJECT_STORE_SUBDIR],
        )
        print(*sorted(pair[0] for pair in it), sep="\n")


def _register_lsversions_command(subparsers):
    parser_ls = subparsers.add_parser("ls-versions", help="list the versions")
    parser_ls.set_defaults(func=lsversions_command)


def getpath_command(args, tool_config, data_store_config):
    if args.data_dir:
        print(data_store_config.data_dir)
    elif args.data_store:
        print(data_store_config.store_location)
    else:
        assert args.vdata is not _UNSPECIFIED  # sanity check!
        if args.vdata is None:
            version = tool_config.grackle_version
        else:
            version = args.vdata
        print(os.path.join(data_store_config.store_location, version))


def _register_getpath_command(subparsers):
    parser_getpath = subparsers.add_parser(
        "getpath",
        description=(
            "Provides the expected filesystem location for data. This command "
            "doesn't care about whether the filesystem location actually exists."
        ),
        help="show expected filesystem location for data.",
    )
    getpath_spec_grp = parser_getpath.add_argument_group(
        title="Target",
        description="specifies the target that we retrieve the path for.",
    ).add_mutually_exclusive_group(required=True)
    getpath_spec_grp.add_argument(
        "--data-dir", action="store_true", help="get path to the data directory"
    )
    getpath_spec_grp.add_argument(
        "--data-store",
        action="store_true",
        help="get path to the data-store (for the protocol version used by this tool)",
    )
    getpath_spec_grp.add_argument(
        "--vdata",
        default=_UNSPECIFIED,
        nargs="?",
        help=(
            "get path to the directory of files-references associated with the "
            "specified version. This command assumes that the version-data was "
            "managed by a version of this tool that uses the same protocol version "
            "as the version returned by --version-protocol. If no version is "
            "specified, it uses the version associated with the --version-dir flag."
        ),
    )
    parser_getpath.set_defaults(func=getpath_command)


def showknownreg_command(args, tool_config, data_store_config):
    f = data_store_config.file_registry_file
    if isinstance(f, io.IOBase):
        lines = f.readlines()
    else:
        with open(data_store_config.file_registry_file, "r") as f:
            lines = f.readlines()
    contents = [
        line for line in lines if len(line.strip()) > 0 and not line.startswith("//")
    ]
    print(*contents, sep="", end="")


def _register_showknownreg_command(subparsers):
    parser_showknownreg = subparsers.add_parser(
        "showknownreg",
        help=(
            "prints the pre-registered file registry expected by the current version "
            "of Grackle"
        ),
    )
    parser_showknownreg.set_defaults(func=showknownreg_command)


def _fmt_registry_lines(fname_cksum_pairs, hash_alg):
    length, suffix = len(fname_cksum_pairs), (",\n", "\n")
    return [
        f'{{"{fname}", "{hash_alg}:{cksum}"}}{suffix[(i+1) == length]}'
        for i, (fname, cksum) in enumerate(sorted(fname_cksum_pairs))
    ]


def calcreg_command(args, tool_config, data_store_config):
    # print the properly file registry information (in the proper format that can be
    # used to configure newer versions of Grackle

    # we use listdir since we are targetting 3.3, but we set things up so that we could
    # use os.scandir
    try:
        it = direntry_iter(args.path, ftype="file", mismatch="eager_err")
    except FileNotFoundError:
        raise ValueError(f"{args.path!r} doesn't specify a directory or file")
    except NotADirectoryError:
        it = [(os.path.basename(args.path), args.path)]

    pairs = [(name, calc_checksum(path, args.hash_name)) for name, path in it]

    with contextlib.ExitStack() as stack:
        if args.output is None:
            file = sys.stdout
        else:
            file = stack.enter_context(open(args.output, "w"))

        '''
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
        '''
        print(*_fmt_registry_lines(pairs, args.hash_name), sep="", end="", file=file)


def _register_calcreg_command(subparsers):
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


def help_command(*args, **kwargs):
    # it might be nice to pipe to a pager (specified by PAGER env variable or

    # here is some logic to strip anchors
    # replace the [[[BEGIN:...]]] & [[[END:...]]] anchors
    _open, _close = r"\[\[\[", r"\]\]\]"
    section_start_anchor = re.compile(
        rf"^{_open}BEGIN-SECTION:([-,_+.! 0-9A-Za-z]+){_close}[ \t]*$"
    )
    generic_anchor = re.compile(rf"^{_open}[-:,_+.! 0-9A-Za-z]+{_close}[ \t]*$")

    for line in _EXTENDED_DESCRIPTION.splitlines():
        m = section_start_anchor.match(line)
        if m:
            section_name = m.group(1)
            print(section_name, len(section_name) * "-", sep="\n")
        elif generic_anchor.match(line):
            continue
        else:
            print(line)


def _register_help_command(subparsers):
    parser_help = subparsers.add_parser(
        "help", help="Display detailed help information about this tool"
    )
    parser_help.set_defaults(func=help_command)


def _add_program_prop_query(parser, flag, value, short_descr):
    """
    add a flag to parser to trigger a control flow that:
    1. shows a fundamental piece of information about the command line program (like a
        version number or the ``--help`` option)
    2. then immediately exits the program
    """

    class _Action(argparse.Action):
        def __call__(self, *args, **kwargs):
            print(value)
            sys.exit(0)

    parser.add_argument(
        flag,
        metavar="",
        action=_Action,
        nargs=0,
        help=f"show associated {short_descr} and exit",
    )


def build_parser(tool_config, prog_name):
    parser = argparse.ArgumentParser(
        prog=prog_name,
        description=(
            "This is a management system for Grackle's data files. Subcommands are "
            "provided to fetch data files, list all available data, and delete data"
        ),
        epilog=f"Invoke `{prog_name} help` to get a detailed overview of the tool",
    )

    # This is a hidden argument. It is only used for sake of testing (we may remove it
    # any time in the future)
    parser.add_argument(
        "--testing-override-registry-file",
        help=argparse.SUPPRESS,  # hides the help message
        default=argparse.SUPPRESS,  # adds no attribute if option wasn't specified
    )
    parser.add_argument(
        "--testing-override-version-grackle",
        help=argparse.SUPPRESS,  # hides the help message
        default=argparse.SUPPRESS,  # adds no attribute if option wasn't specified
    )

    query_l = [
        ("--version-grackle", "Grackle version", tool_config.grackle_version),
        (
            "--version-protocol",
            "data-store protocol version",
            tool_config.protocol_version,
        ),
        ("--cksum-alg", "name of the checksum algorithm", tool_config.checksum_kind),
    ]
    for flag, short_descr, val in query_l:
        _add_program_prop_query(parser, flag, val, short_descr)

    subparsers = parser.add_subparsers(required=True)

    _register_fetch_command(subparsers)
    _register_rm_command(subparsers)
    _register_lsversions_command(subparsers)
    _register_getpath_command(subparsers)
    _register_showknownreg_command(subparsers)
    _register_calcreg_command(subparsers)
    _register_help_command(subparsers)

    return parser


def main(tool_config, data_store_config, prog_name, *, args=None):
    """
    Launch the command

    Returns
    -------
    int
        Specified the exit code
    """
    parser = build_parser(tool_config, prog_name)
    args = parser.parse_args(args=args)

    # handle testing overrides
    if hasattr(args, "testing_override_registry_file"):
        # _replace makes a copy & in the copy any specified attributes are overridden
        data_store_config = data_store_config._replace(
            file_registry_file=args.testing_override_registry_file
        )
    if hasattr(args, "testing_override_version_grackle"):
        # _replace makes a copy & in the copy any specified attributes are overridden
        tool_config = tool_config._replace(
            grackle_version=args.testing_override_version_grackle
        )

    try:
        args.func(args, tool_config=tool_config, data_store_config=data_store_config)
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
        return 78  # https://www.man7.org/linux/man-pages/man3/sysexits.h.3head.html
    except GenericToolError as err:
        print(f"ERROR: {err.args[0]}")
        return 70  # https://www.man7.org/linux/man-pages/man3/sysexits.h.3head.html
    except BaseException:
        print("Unexpected error:", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        return 70  # https://www.man7.org/linux/man-pages/man3/sysexits.h.3head.html
    else:
        return 0


def _default_data_store_config(tool_config, file_registry_file):
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
        default_fetcher=Fetcher.configure_GitHub_url(
            data_repository_url=_REPO_URL,
            contemporaneous_git_hash=_CONTEMPORANEOUS_COMMIT_HASH,
        ),
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
    data_store_config = _default_data_store_config(tool_config, file_registry_file)
    return tool_config, data_store_config


# Here, we define machinery employed when used as a standalone program
# ====================================================================

# to support installing this file as a standalone program, we will need to introduce the
# following procedure to the build-system:
#    - treat this file as a template-file and configure it with CMake's
#      ``configure_file`` command (or invoke ``configure_file.py`` under the classic
#      build system) in order to substitute the names enclosed by the "at sign" symbol
#    - make resulting file executable (and maybe drop the .py suffix)
#    - install it into the bin directory alongside the grackle libraries

if __name__ == "__main__":
    _GRACKLE_VERSION = "@_GRDATA_GRACKLE_VERSION@"
    _FILE_REGISTRY_CONTENTS = """\
@_GRDATA_FILE_REGISTRY_CONTENTS@
"""

    def _check_substitution_problems(var_name, var_value):
        # we use unicode escape sequence, \u0040, that python automatically converts
        # to the "at sign" to prevent the configure_file.py script (used by Grackle's
        # Grackle's build-system) from falsely reporting an error
        if (
            (var_name in var_value)
            or ("\u0040" in var_value)
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
    sys.exit(main(*_CONFIG_PAIR, prog_name="grdata"))
