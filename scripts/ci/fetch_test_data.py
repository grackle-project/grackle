#!/usr/bin/env python3
"""
Download the sample data
"""

# for portability, we restrict ourselves to built-in modules provided with python 3.6+
import importlib.util
import os
import shutil
import sys
import tempfile


def _import_standalone_module(dir_path, module_name):
    path = os.path.join(dir_path, f"{module_name}.py")
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module

_LOCAL_DIR = os.path.dirname(__file__)
grdata = _import_standalone_module(
    os.path.join(_LOCAL_DIR, "..", "..", "src", "python", "pygrackle", "utilities"),
    "grdata"
)


def download_file(url, *, dst=None, dst_dir=None, silent=None, cksum=None):
    """download the file from url to dst"""

    basename = os.path.basename(url if dst is None else dst)
    if (dst is not None) and (dst_dir is not None):
        raise ValueError("dst and dst_dir can't both be specified")
    elif dst is None:
        dst = os.path.join("." if dst_dir is None else dst_dir, basename)

    silent = silent if silent is not None else "true" == os.getenv("CI", None)

    if not silent:
        print(f"> downloading {url}")

    grdata._retrieve_url(url=url, dst=dst, silent=silent)
    if cksum is not None:
        alg_name, _ = cksum.split(":")
        actual = grdata.calc_checksum(fname=dst, alg_name=alg_name)
        if cksum != actual:
            raise RuntimeError(
                f"cksum mismatch for {basename}\n  expected={cksum}\n  actual={actual}"
            )
    return dst


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as tmpdirname:
        envvar = "YT_DATA_DIR"
        data_dir = os.getenv(envvar, None)
        if data_dir is None:
            raise RuntimeError(f"`{envvar}` env var isn't defined")
        elif not os.path.isabs(data_dir):
            raise RuntimeError(f"`{envvar}` env var isn't an absolute path")
        elif not os.path.isdir(data_dir):
            os.mkdir(data_dir)

        url = "https://yt-project.org/data/IsolatedGalaxy.tar.gz"
        ck = "sha256:fc081bd4420efd02f7ba2db7eaf4bff0299d5cc4da93436be55a30c219daaf18"
        archive_path = download_file(url=url, dst_dir=tmpdirname, cksum=ck)
        shutil.unpack_archive(archive_path, data_dir)
    sys.exit(0)
