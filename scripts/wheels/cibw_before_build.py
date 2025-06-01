# the goal is for this script to be very portable

import argparse
import hashlib
import os
import platform
import shutil
import subprocess
import sys

_LOCAL_DIR = os.path.dirname(__file__)
_IS_MACOS = platform.system() == "Darwin"

def _run(*args, check=True, **kwargs):
    print(">", *args, sep=" ", flush=True)
    return subprocess.run(args, check=check, **kwargs)


def calc_checksum(fname, *, alg_name, chunksize=8192):
    """Calculate the checksum for a given fname"""

    # logic is duplicated by scripts in PR #235 AND #307

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
    return f"{alg_name.lower()}:{hash_obj.hexdigest()}"


def download_file(url, *, dst=None, dst_dir=None, quiet=None, cksum=None):
    """download the file from url to dst"""

    # logic is duplicated by scripts in PR #235 AND #307

    basename = os.path.basename(url if dst is None else dst)
    if (dst is not None) and (dst_dir is not None):
        raise ValueError("dst and dst_dir can't both be specified")
    elif dst is None:
        dst = os.path.join("." if dst_dir is None else dst_dir, basename)

    quiet = quiet if quiet is not None else "true" == os.getenv("CI", None)

    choices = [
        ("wget", None, "--output-document", "--quiet"),
        ("curl", "-L", "--output", "--silent"),
    ]
    for tool, loc_flag, outflag, quietflag in choices:
        if shutil.which(tool) is not None:
            quiet_args = [quietflag] if quiet else []
            url_args = [url] if loc_flag is None else [loc_flag, url]
            _cmd = [tool, outflag, dst] + quiet_args + url_args
            _run(*_cmd, check=True)

            if cksum is not None:
                alg_name, _ = cksum.split(":")
                actual = calc_checksum(fname=dst, alg_name=alg_name)
                if cksum != actual:
                    raise RuntimeError(
                        f"cksum mismatch for {basename}\n"
                        f"  expected={cksum}\n  actual={actual}"
                    )
            return dst
    raise RuntimeError("neither wget nor curl is installed")


def get_gfortran(depend_dir):
    if not _IS_MACOS:
        return None  # gfortran already exists on Linux
    print("downloading gfortran", flush=True)

    release, _, machine = platform.mac_ver()

    # both branches are inspired by scipy's tool/wheels/cibw_before_build_macos.sh
    if machine == "arm64":
        dmg_path = download_file(
            url="https://github.com/fxcoudert/gfortran-for-macOS/releases/download/12.1-monterey/gfortran-ARM-12.1-Monterey.dmg",
            cksum="sha256:e2e32f491303a00092921baebac7ffb7ae98de4ca82ebbe9e6a866dd8501acdf",
            dst=f"{depend_dir}/gfortran.dmg",
        )

        # install the disk image
        _run("hdiutil", "attach", "-mountpoint", "/Volumes/gfortran", dmg_path)
        _run(
            "sudo",
            "installer",
            "-pkg",
            "/Volumes/gfortran/gfortran.pkg",
            "-target",
            "/",
        )
        _run("type", "-p", "gfortran")
    else:
        tgz_path = download_file(
            url="https://github.com/isuruf/gcc/releases/download/gcc-11.3.0-2/gfortran-darwin-x86_64-native.tar.gz",
            cksum="sha256:981367dd0ad4335613e91bbee453d60b6669f5d7e976d18c7bdb7f1966f26ae4",
            dst=f"{depend_dir}/gfortran.tar.gz",
        )

        if not os.path.isdir("/opt"):
            _run("sudo", "mkdir", "/opt/")
        _run("sudo", "tar", "-xv", "-C", "/opt", "-f", tgz_path)

        # link the following libraries into /usr/local so we don't need to provide the
        # -L and -rpath flags to give hints to the static linker and dynamic loader
        libs = [
            "libgfortran.dylib",
            "libgfortran.5.dylib",
            "libgcc_s.1.dylib",
            "libgcc_s.1.1.dylib",
            "libquadmath.dylib",
            "libquadmath.0.dylib",
        ]
        for lib in libs:
            os.symlink(
                f"/opt/gfortran-darwin-x86_64-native/lib/{lib}", f"/usr/local/lib/{lib}"
            )
        os.symlink(
            "/opt/gfortran-darwin-x86_64-native/bin/gfortran", "/usr/local/bin/gfortran"
        )


def get_hdf5(depend_dir, compile_hl_api=False, build_type="Release"):
    # this logic could ostensibly be placed into src/python/CMakeLists.txt (to be
    # invoked by scikit-build-core)
    # -> if we do that, we need to make it **VERY** clear not to use the logic within
    #    the CMake build of the core grackle-library (people almost never want that
    #    because it means that downstream simulation codes can't link against hdf5)
    # -> The advantage to keeping this separate is we can avoid recompiling hdf5 each
    #    time we compile a wheel with a different ABI (py310, py311) on a given platform
    # -> the practice of distributing the optional libz dependency with the wheel on a
    #    subset of platforms makes complicates transcription of this logic into CMake
    # -> we would also need to make sure that hdf5 and zlib (where applicable) are
    #    properly renamed & put into the appropriate directories when "repairing" the
    #    wheel with auditwheel/delocate

    # we decided that we want to ship hdf5 with builtin gzip compression support
    # (like h5py, we don't ship SZIP support due to patent considerations)

    # 1. Pre-download & unpack the archive for HDF5
    h5_archive_path = download_file(
        url="https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5_1.14.6.tar.gz",
        cksum="sha256:09ee1c671a87401a5201c06106650f62badeea5a3b3941e9b1e2e1e08317357f",
        dst_dir=depend_dir,
    )
    h5_srcdir, h5_builddir = [os.path.join(depend_dir, p) for p in ("h5src", "h5build")]
    os.mkdir(h5_srcdir)
    _run("tar", "-xf", h5_archive_path, "-C", h5_srcdir, "--strip-components=1")

    # 2. handle the zlib dependency (when we include it in the precompiled wheel)
    #    - In an earlier version, I spent a little time trying to use the machinery of
    #      the HDF5 build-system to try to compile and build zlib as part of building
    #      HDF5. I was trying to manually download the file, so we could validate the
    #      checksum (the HDF5 can't currently do that!), but I couldn't get it to work
    #    - thus, we manually install it ourselves
    if sys.platform.startswith("win32"):
        raise RuntimeError("we need to implement logic to download & install zlib")
    elif _IS_MACOS:
        zlib_version, zlib_url, zlib_cksum = (
            "1.3.1",
            "https://github.com/madler/zlib/releases/download/v1.3.1/zlib-1.3.1.tar.gz",
            "sha256:9a93b2b7dfdac77ceba5a558a580e74667dd6fede4585b91eefb60f03b72df23"
        )
        assert f"/v{zlib_version}/" in zlib_url # sanity-check!
        with open(f"{h5_srcdir}/config/cmake/ZLIB/CMakeLists.txt", "r") as f:
            assert f'set(VERSION "{zlib_version}")' in f.read() # sanity-check!
        zlib_archive_path = download_file(
            url=zlib_url, cksum=zlib_cksum, dst_dir=depend_dir
        )
        # create the list of zlib-related args to pass to CMake
        #external_zlib_args = [
        #    "-DHDF5_ALLOW_EXTERNAL_SUPPORT=TGZ",
        #    f"-DZLIB_TGZ_NAME={os.path.abspath(zlib_archive_path)}",
        #    f"-DTGZPATH={os.path.abspath(zlib_archive_path)}"
        #]

        external_zlib_args = [
            "-DHDF5_ALLOW_EXTERNAL_SUPPORT=TGZ",
            f"-DZLIB_TGZ_NAME={os.path.basename(zlib_archive_path)}",
            f"-DZLIB_TGZ_ORIGPATH={os.path.abspath(zlib_archive_path)}",
            f"-DZLIB_USE_LOCALCONTENT=ON"
        ]
    else:
        # in this case, we use zlib installation provided by the
        external_zlib_args = []

    # 3. Configure, Compile, and Install HDF5
    _run(  # configure the build
        "cmake",
        f"-S{h5_srcdir}",
        f"-B{h5_builddir}",
        f"-DCMAKE_BUILD_TYPE={build_type}",
        "-DCMAKE_INSTALL_PREFIX=/usr/local/",
        "-DBUILD_STATIC_LIBS=OFF",
        # we only compile the high-level api if we want to be able to compile h5py on
        # platforms (e.g. musllinux) where h5py doesn't provide binary wheels
        f"-DHDF5_BUILD_HL_LIB={compile_hl_api!s}",
        # here we just disable a bunch of unneeded components
        "-DBUILD_TESTING=OFF",
        "-DHDF5_BUILD_CPP_LIB=OFF",
        "-DHDF5_BUILD_EXAMPLES=OFF",
        "-DHDF5_BUILD_FORTRAN=OFF",
        "-DHDF5_BUILD_JAVA=OFF",
        "-DHDF5_BUILD_TOOLS=OFF",
        "-DHDF5_BUILD_PARALLEL_TOOLS=OFF",
        "-DHDF5_BUILD_STATIC_TOOLS=OFF",
        "-DHDF5_ENABLE_SZIP_SUPPORT=OFF",
        "-DHDF5_USE_ZLIB_STATIC=OFF",
        "-DHDF5_ENABLE_Z_LIB_SUPPORT=ON",  # the option-name changes in HDF5 2.0
        *external_zlib_args
    )
    _run("cmake", "--build", h5_builddir)  # compile hdf5
    if _IS_MACOS:
        _run("sudo", "cmake", "--install", h5_builddir)
    else:
        _run("cmake", "--install", h5_builddir)


def handle_license(project_dir):
    prefix = "pygrackle/.dylibs" if _IS_MACOS else "pygrackle.libs"
    template_path = os.path.join(_LOCAL_DIR, "binary_license_annex.txt.in")
    print(f"replacing @PREFIX@ with {prefix} in {template_path}")
    configure_bin = os.path.join(_LOCAL_DIR, "..", "configure_file.py")
    annex_path = os.path.join(_LOCAL_DIR, "_binary_license_annex.txt")
    literal_linenos = ["112", "113", "902"]
    _run(
        configure_bin,
        f"--input={template_path}",
        f"--output={annex_path}",
        "--clobber",
        f"PREFIX={prefix}",
        "--literal-linenos",
        *literal_linenos,
    )

    annex_list = [annex_path]
    if _IS_MACOS:
        annex_list.append(os.path.join(_LOCAL_DIR, "zlib_license_for_macos.txt"))

    license_path = os.path.join(project_dir, "LICENSE")
    assert os.path.isfile(license_path)
    print(f"append binary-license-annex to {license_path}")
    with open(license_path, "a") as fout:
        for path in annex_list:
            with open(path, "r") as fsrc:
                shutil.copyfileobj(fsrc, fout)


parser = argparse.ArgumentParser()
parser.add_argument(
    "--compile-hl-h5", action="store_true", help="compile high-level h5 api"
)
parser.add_argument("depend_dir", help="scratch-space while setting up dependencies")
parser.add_argument("project_dir", help="path to the project-dir")

if __name__ == "__main__":
    if sys.platform.startswith("win32"):
        raise ValueError("there's currently no support for windows")
    args = parser.parse_args()
    depend_dir, project_dir = args.depend_dir, args.project_dir

    print("Showing Environment")
    for key, value in os.environ.items():
        print(f"  {key!s}={value!s}")
    print("", end="\n", flush=True)

    handle_license(args.project_dir)

    # make the dependency-directory
    print(f"creating: {depend_dir}",flush=True)
    if os.path.isdir(args.depend_dir):
        print("the directory already exists -- nothing needs to be done", flush=True)
    else:
        os.mkdir(args.depend_dir)
        # install gfortran (does nothing on linux, where gfortran is already installed)
        get_gfortran(args.depend_dir)
        # download, build, and install HDF5
        get_hdf5(args.depend_dir, compile_hl_api=args.compile_hl_h5)
