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
import sys
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
from common import (
    LinuxDistroInfo,
    ScriptError,
    configure_logger,
    exec_cmd,
    logger,
    export_env_var_assignment,
)


class Conf(NamedTuple):
    manual_install_dir: str
    dry_run: bool
    os_str: str
    distro_info: Optional[LinuxDistroInfo]
    mutate_path: bool = True

    @classmethod
    def from_args(cls, args: argparse.Namespace):
        os_str = platform.system().lower()
        distro_info = LinuxDistroInfo.infer_from_system() if os_str == "linux" else None
        return Conf(
            manual_install_dir=os.path.expanduser(args.manual_install_dir),
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
    dst_basename: str


class APTPackage(NamedTuple):
    name: str


class CodeSource(NamedTuple):
    """Represents the source for a piece of software"""

    source: Union[APTPackage, Callable[[], None]]
    setup_fn: Optional[Callable[[Conf], None]] = None
    symlinks: Optional[List[SymlinkSpec]] = None


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


def _extract_major_minor(version: str) -> Tuple[int, int]:
    m = re.match(r"^(0|[1-9]\d*)\.(0|[1-9]\d*)(?:[^0-9].*)?$", version)
    if m is None:
        raise ScriptError(f"{version} doesn't have a major and minor version")
    return tuple(int(e) for e in m.groups())


def _download_and_install_precompiled(
    url_prefix: str, cksum_fname: str, desired_fname: str, outdir: str, cksum_kind: str
):
    # determine the checksum (to let us check for errors in the download)
    # -> this doesn't really add much security (if someone compromises the
    #    precompiled binary, they probably also compromised the cksum file)
    logger.info(f"fetching the checksum for {desired_fname}")
    with tempfile.TemporaryFile() as fp:
        download_file(url=f"{url_prefix}/{cksum_fname}", dst=fp)
        fp.seek(0)
        cksum = _search_cksum_file(fp, target_fname=desired_fname)
    if cksum is None:
        raise RuntimeError(f"could not find the checksum for {desired_fname}")
    elif ":" not in cksum:
        cksum = f"{cksum_kind}:{cksum}"

    # download the full file
    logger.info(f"downloading {desired_fname}")
    with tempfile.NamedTemporaryFile() as fp:
        tmp_path = fp.name
        download_file(url=f"{url_prefix}/{desired_fname}", dst=tmp_path, cksum=cksum)
        logger.info("unpacking the contents")
        exec_cmd("tar", "-x", "--strip-components", "1", "-C", outdir, "-f", tmp_path)


def get_cmake_source(conf: Conf, version: Optional[str]) -> CodeSource:
    if version is None:
        return CodeSource(source=APTPackage("cmake"))
    else:
        # in this case, we try to download a precompiled binary
        # (building from source is also very easy, but will be slower)
        try:
            major, minor = _extract_major_minor(version)
        except ScriptError:
            raise ScriptError(f"{version} isn't a valid version for CMake") from None
        url_prefix = f"https://cmake.org/files/v{major}.{minor}"

        if (major, minor) < (3, 20):
            raise RuntimeError(f"CMake url pattern isn't handled for {major}.{minor}")
        elif f"{major}.{minor}" == version:
            logger.info(f"inferring full CMake version from {version}")
            with tempfile.TemporaryFile() as fp:
                download_file(url=f"{url_prefix}/cmake-latest-files-v1.json", dst=fp)
                fp.seek(0)
                version = json.load(fp)["version"]["string"]

        outdir = os.path.join(conf.manual_install_dir, f"cmake-{version}")

        def install_cmake():
            cksum_fname = f"cmake-{version}-SHA-256.txt"
            if conf.os_str == "darwin":
                desired = f"cmake-{version}-macos-universal.dmg"
            else:
                desired = f"cmake-{version}-{conf.os_str}-{platform.machine()}.tar.gz"
            if not conf.dry_run:
                if not os.path.isdir(outdir):
                    os.mkdir(outdir)
                _download_and_install_precompiled(
                    url_prefix, cksum_fname, desired, outdir, cksum_kind="sha256"
                )

        installed_bin = os.path.join(outdir, "bin", "cmake")
        return CodeSource(
            source=install_cmake,
            symlinks=[
                SymlinkSpec(src=installed_bin, dst_basename=f"cmake-{version}"),
                SymlinkSpec(src=installed_bin, dst_basename="cmake"),
            ],
        )


def setup_llvm_repository(conf: Conf, llvm_version: int):
    """This function registers LLVM's official apt repository"""
    # this logic is loosely based on the from https://apt.llvm.org/llvm.sh

    _descr = "to register LLVM's official apt repository"
    if conf.distro_info is None:
        raise ScriptError(f"{_descr}: linux distribution info must be known")
    elif conf.distro_info.id not in ["ubuntu", "debian"]:
        raise ScriptError(f"{_descr}: only possible on debian/ubuntu")
    elif conf.distro_info.codename is None:
        raise ScriptError(f"{_descr}: requires knowledge of the codename")

    logger.info("Fetch and register key for authenticating LLVM's APT repository")
    gpg_dst = "/etc/apt/trusted.gpg.d/apt.llvm.org.asc"
    gpg_src_url = "https://apt.llvm.org/llvm-snapshot.gpg.key"
    with tempfile.NamedTemporaryFile() as fp:
        tmp_path = fp.name
        download_file(gpg_src_url, dst=tmp_path)
        exec_cmd("sudo", "cp", tmp_path, gpg_dst, dry_run=conf.dry_run)
        exec_cmd("sudo", "chmod", "a+r", gpg_dst, dry_run=conf.dry_run)

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


def get_clang_tidy_source(conf: Conf, version: Optional[str]) -> CodeSource:
    if version is None:
        raise ScriptError("you must associate clang-tidy with a version")
    return CodeSource(
        source=APTPackage(f"clang-tidy-{version}"),
        setup_fn=functools.partial(
            setup_llvm_repository, llvm_version=version, conf=conf
        ),
        symlinks=[
            SymlinkSpec(src=f"/usr/bin/clang-tidy-{version}", dst_basename="clang-tidy")
        ],
    )


class SimpleSourceFactory(NamedTuple):
    """Used when APT (with the default repository) is the only way install software"""

    name: str
    apt_pkg_name: str

    def __call__(self, conf: Conf, version: Optional[str]) -> CodeSource:
        if version is not None:
            raise ScriptError(f"{self.name!s} cannot be specified with a version")
        return CodeSource(source=APTPackage(name=self.apt_pkg_name))


_SOFTWARE: Dict[str, Callable[[Conf, str], CodeSource]] = {
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
    preset_set = set(presets)  # <- avoid mutating presets
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


def _apt_install(conf: Conf, names: Collection[str]):
    assert len(names) > 0, "sanity check failed"
    exec_cmd("sudo", "apt-get", "update", dry_run=conf.dry_run)
    exec_cmd("sudo", "apt-get", "install", "-y", *sorted(names), dry_run=conf.dry_run)


def _handle_symlinks(conf: Conf, symlink_sets: List[Tuple[str, List[SymlinkSpec]]]):
    bin_dir = os.path.abspath(os.path.join(conf.manual_install_dir, "bin"))
    logger.info(f"  ensure BINDIR={bin_dir} exists")
    if not conf.dry_run:
        os.makedirs(bin_dir, exist_ok=True)
        export_env_var_assignment(f"export PATH={bin_dir}:$PATH")
    for descr, symlink_spec_set in symlink_sets:
        for spec in symlink_spec_set:
            msg = f"  for {descr}: SRC={spec.src}, DST=$BINDIR/{spec.dst_basename}"
            logger.info(msg)
            if conf.dry_run:
                continue
            elif not os.path.isfile(spec.src):
                raise RuntimeError(
                    f"can't make a symlink to {spec.src}, (it doesn't seem to exist)"
                )
            else:
                os.symlink(spec.src, os.path.join(bin_dir, spec.dst_basename))


def install_software(software_str_itr: Iterable[str], conf: Conf):
    """Does all of the heavy-lifting"""
    # software strings may occur more than once in software_str_itr

    logger.info(f"Sys Info: {conf.sys_info_str()}")

    if not conf.dry_run:
        try:
            os.mkdir(conf.manual_install_dir)
        except FileExistsError:
            pass  # <- this case is ok!
        except OSError:
            raise ScriptError(
                f"{conf.manual_install_dir} doesn't exist & can't easily be created"
            ) from None

    # fill all_software with name : (version, software_name) entries
    all_software = {}
    for software_str in software_str_itr:
        name, new_v = _get_name_version(software_str)
        if name not in _SOFTWARE:
            raise ScriptError(f"{software_str} doesn't specify known software")
        old_v, _ = all_software.get(name, (None, None))
        if (old_v is not None) and (new_v is None):
            continue
        elif (old_v is not None) and (old_v != new_v):
            logger.warning(f"override {name!r} version: use {new_v} instead of {old_v}")
        all_software[name] = (new_v, software_str)

    # infer the required steps
    cmd_list: List[Tuple[str, str, Callable[[], None]]] = []
    apt_packages: Set[str] = set()
    symlink_sets: List[Tuple[str, List[SymlinkSpec]]] = []

    for name, (version, software_str) in all_software.items():
        src = _SOFTWARE[name](conf=conf, version=version)
        if src.setup_fn is not None:
            cmd_list.append(("setup step", software_str, src.setup_fn))
        if callable(src.source):
            cmd_list.append(("manual step", software_str, src.source))
        else:
            assert isinstance(src.source, APTPackage), "sanity check"
            apt_packages.add(src.source.name)
        if src.symlinks is not None:
            symlink_sets.append((software_str, src.symlinks))

    if len(apt_packages) > 0:
        fn = functools.partial(_apt_install, conf=conf, names=apt_packages)
        cmd_list.append(("invoke APT", None, fn))

    if len(symlink_sets) > 0:
        fn = functools.partial(_handle_symlinks, conf=conf, symlink_sets=symlink_sets)
        cmd_list.append(("Handle Symlinks (And Paths)", None, fn))

    # execute the commands
    n_steps = len(cmd_list)
    for i, (descr, target, cmd) in enumerate(cmd_list, start=1):
        msg = f"{descr} for {target}" if target is not None else descr
        logger.info(f"Dependency Install {i}/{n_steps}: {msg}")
        cmd()


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
    sys.exit(main(parser.parse_args()))
