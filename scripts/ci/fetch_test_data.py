#!/usr/bin/env python3
"""
Download the sample data
"""

# for portability, we restrict ourselves to built-in modules provided with python 3.6+
import contextlib
import hashlib
from math import log10
import os
import shutil
import sys
import tempfile
import urllib.request


_CHUNKSIZE = 8192  # default chunksize used for file operations


# duplicates logic in PR #235
def _download_status_updater(tot_bytes, *, silent=False):
    """provides information about the current download"""
    power_div_3 = int(log10(tot_bytes) // 3) if tot_bytes > 0 else 0
    div, unit = 1000.0**power_div_3, (" B", "KB", "MB", "GB")[power_div_3]
    done = tot_bytes <= 0 or silent or not sys.stdout.isatty()

    def _update(size):
        nonlocal done
        m = "\n" if size is None else f"\r{size / div:.2f}/{tot_bytes / div:.2f} {unit}"
        if not done:
            print(m, end="", flush=True)
            done = size is None

    return _update


# duplicates logic in PR #235
def _retrieve_url(url, dst, *, silent=False, chunksize=_CHUNKSIZE):
    """Download the datafile from url to dst."""
    req = urllib.request.Request(url)
    with contextlib.ExitStack() as stack:
        out_file = stack.enter_context(open(dst, "wb"))
        response = stack.enter_context(urllib.request.urlopen(req))
        total_bytes = int(response.headers.get("Content-Length", -1))
        update_progress = _download_status_updater(total_bytes, silent=silent)
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


# duplicates logic in PR #235
def calc_checksum(fname, *, alg_name, chunksize=_CHUNKSIZE):
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
    return f"{alg_name.lower()}:{hash_obj.hexdigest()}"


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

    _retrieve_url(url=url, dst=dst, silent=silent)
    if cksum is not None:
        alg_name, _ = cksum.split(":")
        actual = calc_checksum(fname=dst, alg_name=alg_name)
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
