#!/usr/bin/env python3
"""
This program seeks to install continuous integration dependencies.

At the moment, this should work for ubuntu/debian images.

Longer term, the goals are to:
- unify this logic with the cibw scripts (that means we need to add support for macOS)
- make it possible to install arbitrary versions of different compiler toolchains
  and different versions of dependencies. For example, it would be good practice
  for us to perform a test build with the oldest compatible versions of
  hdf5/cmake/etc.
"""

# for portability, we restrict ourselves to built-in modules provided with python 3.6

import argparse
import functools
import itertools
import json
import os
import platform
import re
import tempfile
from typing import (
    Dict,
    Callable,
    Collection,
    IO,
    Iterable,
    List,
    NamedTuple,
    Optional,
    Set,
    Tuple,
    Union,
)

from fetch_test_data import download_file
from common import LinuxDistroInfo, ScriptError, configure_logger, exec_cmd, logger


class Conf(NamedTuple):
    manual_install_dir: str
    dry_run: bool
    os_str: str
    distro_info: Optional[LinuxDistroInfo]

    @classmethod
    def from_args(cls, args: argparse.Namespace):
        os_str = platform.system().lower()
        distro_info = LinuxDistroInfo.infer_from_system() if os_str == "linux" else None
        return Conf(
            manual_install_dir=args.manual_install_dir,
            dry_run=args.dry_run,
            os_str=os_str,
            distro_info=distro_info,
        )

    def sys_info_str(self) -> str:
        tmp = f"os={self.os_str}"
        distro_info = self.distro_info._asdict() if self.distro_info else {}
        return tmp + "".join(f", {k}={v}" for k, v in distro_info.items())


class SymlinkSpec(NamedTuple):
    src: str
    dst: str


class CodeSourceAPT(NamedTuple):
    """Encapsulates information for installing a package via APT"""

    pkg_name: str
    # the following attributes are usually just used with non-default APT repositories
    setup_fn: Optional[Callable[[Conf], None]] = None
    symlinks: Optional[List[SymlinkSpec]] = None


# The source of each piece of software is associated with either
#   1. an APT package OR
#   2. a simple callable function that manually performs the installation itself
# Given a version string (or `None`) a source factory returns one of the above 2 choices
SimpleCallable = Callable[[Conf], Optional[List[SymlinkSpec]]]
SoftwareSource = Union[CodeSourceAPT, SimpleCallable]
SourceFactory = Callable[[Optional[str]], SoftwareSource]


def _search_cksum_file(fp: IO[str], target_fname: str) -> Optional[str]:
    """
    Extracts the checksum from a file using the standard format (or returns
    None if the target_filename can't be found
    """
    # for context, each line of such a file holds `<cksum>  <filename>`
    for line in fp:
        cksum, fname = line.rstrip().decode("utf-8").split()
        if fname == target_fname:
            return cksum


def install_cmake(target_version: str, conf: Conf) -> List[SymlinkSpec]:
    """Download an install a target cmake version"""
    # It's also easy to build from source, but this is faster

    # extract the prefix of target_version that holds major and minor version numbers
    m = re.match(r"^(\d+\.\d+).*", target_version)
    if m is None:
        raise ScriptError(f"{target_version} is an invalid CMake version")
    major_minor_version = m.group(1)
    major, minor = [int(e) for e in major_minor_version.split(".")]
    if (major, minor) < (3, 20):
        raise RuntimeError(f"the url pattern isn't handled for {major_minor_version}")

    # set up some basic properties
    prefix = "https://cmake.org/files"
    url_dir = f"{prefix}/v{major_minor_version}"
    outdir = os.path.join(conf.manual_install_dir, f"cmake-{major_minor_version}")

    # if the caller specified "major.minor", let's infer "major.minor.patch"
    if major_minor_version == target_version:
        # we are missing a patch version, infer the latest version
        logger.info(f"inferring newest CMake patch associated with {target_version}")
        with tempfile.TemporaryFile() as fp:
            download_file(url=f"{url_dir}/cmake-latest-files-v1.json", dst=fp)
            fp.seek(0)
            target_version = json.load(fp)["version"]["string"]

    logger.info(f"prepare to download CMake version: {target_version}")

    # determine the name of the file that will be downloaded
    if conf.os_str == "darwin":
        desired = f"cmake-{target_version}-macos-universal.dmg"
    else:
        arch = platform.machine()
        arch = "x86_64"
        desired = f"cmake-{target_version}-{conf.os_str}-{arch}.tar.gz"

    # determine the sha256 checksum
    # -> this lets us check for errors in the download
    # -> this doesn't really add much security (if someone compromises the
    #    precompiled binary, they probably also compromised the cksum file)
    logger.info(f"fetching the checksum for {desired}")
    with tempfile.TemporaryFile() as fp:
        download_file(url=f"{url_dir}/cmake-{target_version}-SHA-256.txt", dst=fp)
        fp.seek(0)
        cksum = _search_cksum_file(fp, target_fname=desired)
    if cksum is None:
        raise RuntimeError(f"could not find the checksum for {desired}")
    elif ":" not in cksum:
        cksum = "sha256:" + cksum

    # download the full file
    logger.info(f"downloading {desired}")
    with tempfile.NamedTemporaryFile() as fp:
        tmp_path = fp.name
        download_file(url=f"{url_dir}/{desired}", dst=tmp_path, cksum=cksum)
        logger.info("unpacking the contents")
        if not conf.dry_run:
            if not os.path.isdir(outdir):
                os.mkdir(outdir)
            exec_cmd("tar", "-x", "-C", outdir, "-f", tmp_path)

    installed_bin = os.path.join(outdir, "cmake")

    return [
        SymlinkSpec(installed_bin, f"/usr/bin/cmake-{target_version}"),
        SymlinkSpec(installed_bin, f"/usr/bin/cmake-{major_minor_version}"),
        # this last choice is VERY hacky
        # -> we probably want to something clever with paths instead
        SymlinkSpec(installed_bin, "/usr/bin/cmake"),
    ]


def get_cmake_source(version: Optional[str]) -> Union[CodeSourceAPT, SimpleCallable]:
    if version is None:
        return CodeSourceAPT("cmake")
    return functools.partial(install_cmake, version)


def setup_llvm_repository(conf: Conf, llvm_version: int):
    """This function registers LLVM's official apt repository"""
    # this logic is loosely based on the from https://apt.llvm.org/llvm.sh

    _descr = "to register LLVM's official apt repository"
    if conf.distro_info is None:
        raise ScriptError("linux distribution info must be known")
    elif conf.distro_info.id not in ["ubuntu", "debian"]:
        raise ScriptError(f"{_descr}: only possible on debian/ubuntu")
    elif conf.distro_info.codename is None:
        raise ScriptError(f"{_descr}: requires knowledge of the codename")

    logger.info("Fetch and register key for authenticating LLVM's APT repository")
    gpg_dst = "/etc/apt/trusted.gpg.d/apt.llvm.org.asc"
    gpg_src_url = "https://apt.llvm.org/llvm-snapshot.gpg.key"
    if not conf.dry_run:
        download_file(gpg_src_url, dst=gpg_dst)

    # Part 2: record the actual source of the package
    # -> the formula for determining the uri and suite-name is fairly
    #    consistent unless you are describing Debian Testing or the
    #    LLVM version under active development

    logger.info("Register the appropriate LLVM repository with APT")
    uri = f"http://apt.llvm.org/{conf.distro_info.codename}/"
    suite = f"llvm-toolchain-{conf.distro_info.codename}-{llvm_version}"
    cmd = (  # the last arg is intentionally a single string
        "add-apt-repository",
        "--yes",
        "--sourceslist",
        f"deb {uri} {suite} main",
    )
    exec_cmd("sudo", *cmd, dry_run=conf.dry_run)


def get_clang_tidy_source(version: Optional[str]) -> CodeSourceAPT:
    if version is None:
        raise ScriptError("you must associate clang-tidy with a version")
    # the choice to make this symlink is VERY hacky
    tmp = SymlinkSpec(src=f"/usr/bin/clang-tidy-{version}", dst="/usr/bin/clang-tidy")
    return CodeSourceAPT(
        pkg_name=f"clang-tidy-{version}",
        setup_fn=functools.partial(setup_llvm_repository, llvm_version=version),
        symlinks=[tmp],
    )


class SimpleSourceFactory(NamedTuple):
    """Used when APT (with the default repository) is the only way install software"""

    name: str
    apt_pkg_name: str

    def __call__(self, version: Optional[str]) -> CodeSourceAPT:
        if version is not None:
            raise ScriptError(f"{self.name!s} cannot be specified with a version")
        return CodeSourceAPT(pkg_name=self.apt_pkg_name)


_SOFTWARE: Dict[str, SourceFactory] = {
    "castxml": SimpleSourceFactory(name="castxml", apt_pkg_name="castxml"),
    "clang-tidy": get_clang_tidy_source,
    "cmake": get_cmake_source,
    "hdf5": SimpleSourceFactory(name="hdf5", apt_pkg_name="libhdf5-serial-dev"),
    "gfortran": SimpleSourceFactory(name="gfortran", apt_pkg_name="gfortran"),
    "gcc": SimpleSourceFactory(name="gcc", apt_pkg_name="gcc"),
    "g++": SimpleSourceFactory(name="g++", apt_pkg_name="g++"),
    "gmake": SimpleSourceFactory(name="gmake", apt_pkg_name="make"),  # <- GNU Make
    "libtool": SimpleSourceFactory(name="libtool", apt_pkg_name="libtool-bin"),
    "ninja": SimpleSourceFactory(name="ninja", apt_pkg_name="ninja-build"),
}


# the values are a tuple([<dependent_presets>], [<software>])
_PRESETS = {
    "common": ([], ["hdf5", "gfortran", "gcc", "g++"]),
    "cmake-build": (["common"], ["cmake", "ninja"]),
    "classic-build": (["common"], ["gmake", "libtool"]),
    "corelib-tests": (["cmake-build"], ["castxml"]),
    "static-analysis": (["cmake-build"], ["clang-tidy=20"]),
}
_PRESETS["all-presets"] = (list(_PRESETS), [])


def software_from_presets(presets: Iterable[str]) -> Iterable[str]:
    preset_set = set(presets)  # <- avoid mutating args.preset
    while preset_set:
        preset_l, software_l = _PRESETS[preset_set.pop()]
        preset_set.update(preset_l)
        yield from software_l


def _get_name_version(software_str: str) -> Tuple[str, Optional[str]]:
    count = software_str.count("=")
    if count == 0:
        return software_str, None
    elif count == 1:
        return software_str.split("=")
    raise ScriptError(f"{software_str!r} is an invalid string for specifying software")


def apt_install(conf: Conf, names: Collection[str]):
    assert len(names) > 0, "sanity check failed"
    exec_cmd("sudo", "apt-get", "update", dry_run=conf.dry_run)
    exec_cmd("sudo", "apt-get", "install", "-y", *sorted(names), dry_run=conf.dry_run)


def install_software(software_str_itr: Iterable[str], conf: Conf):
    """Does all of the heavy-lifting"""
    # software strings may occur more than once in software_str_itr

    logger.info(f"Sys Info: {conf.sys_info_str()}")

    # fill all_software with name : (version, software_name) entries
    all_software = {}
    for software_str in software_str_itr:
        name, new_v = _get_name_version(software_str)
        if name not in _SOFTWARE:
            raise ScriptError(f"{software_str} doesn't specify known software")
        old_v, _ = all_software.get(name, (None, None))
        if old_v is not None:
            logger.warning(f"override {name!r} version: use {new_v} instead of {old_v}")
        all_software[name] = (new_v, software_str)

    # infer the required steps
    cmd_list: List[Tuple[str, str, SimpleCallable]] = []
    apt_packages: Set[str] = set()
    symlink_sets: List[Tuple[str, List[SymlinkSpec]]] = []

    for name, (version, software_str) in all_software.items():
        src = _SOFTWARE[name](version=version)
        if callable(src):
            cmd_list.append(("manual step", software_str, src))
        else:
            assert isinstance(src, CodeSourceAPT)
            apt_packages.add(src.pkg_name)
            if src.setup_fn is not None:
                cmd_list.append(("setup step", software_str, src.setup_fn))
            if src.symlinks is not None:
                symlink_sets.append((software_str, src.symlinks))

    if len(apt_packages) > 0:
        fn = functools.partial(apt_install, names=apt_packages)
        cmd_list.append(("apt invocations", None, fn))

    n_steps = len(cmd_list) + 1  # add 1 for symlinks

    # actually do the installation
    for i, (descr, target, cmd) in enumerate(cmd_list, start=1):
        msg = f"{descr} for {target}" if target is not None else descr
        logger.info(f"Dependency Install {i}/{n_steps}: {msg}")
        extra_symlinks = cmd(conf)
        if extra_symlinks:
            assert target is not None, "sanity check"
            symlink_sets.append((target, extra_symlinks))

    logger.info(f"Dependency Install {n_steps}/{n_steps}: Make Symlinks (if needed)")
    for descr, symlink_set in symlink_sets:
        for symlink_spec in symlink_set:
            msg = f"  for {descr}: SRC={symlink_spec.src}, DST={symlink_spec.dst}"
            logger.info(msg)
            if not conf.dry_run:
                os.symlink(symlink_spec.src, symlink_spec.dst)


def main(args: argparse.Namespace):
    configure_logger(color=args.color)

    try:
        itr_a = [] if args.preset is None else software_from_presets(args.preset)
        itr_b = [] if args.software is None else args.software
        install_software(itertools.chain(itr_a, itr_b), Conf.from_args(args))

    except ScriptError as err:
        # in this case, we handle "expected errors"
        logger.error(f"ERROR: {err.args[0]}")
        return 70  # https://www.man7.org/linux/man-pages/man3/sysexits.h.3head.html
    except BaseException:
        # here we handle all other exceptions (e.g. programming errors,
        # KeyboardInterrupt). Generally we want a standard traceback in these cases
        logger.error("Unexpected error:")
        raise
    else:
        return 0


_SHORT_DESCR = "Install dependencies for use in continuous integration."
parser = argparse.ArgumentParser(description=_SHORT_DESCR, allow_abbrev=False)
parser.add_argument("--color", action="store_true", help="use color")
parser.add_argument("--dry-run", action="store_true")
parser.add_argument(
    "--software",
    action="append",
    help="""
        Each occurrence adds to the list of software that this tool installs. To specify
        a particular version, write the value as `<name>=<version>`
    """,
)
parser.add_argument(
    "--preset",
    action="append",
    choices=list(_PRESETS),
    help="""
        Each occurrence of this flag appends the associated value to the list of
        software presets that this tool will install
    """,
)
parser.add_argument(
    "--manual-install-dir",
    action="store",
    default="/opt",
    help="specifies destination for manually installed dependencies",
)

if __name__ == "__main__":
    main(parser.parse_args())
