#!/usr/bin/env python3
# download the sample data

import hashlib
import os
import subprocess
import shutil
import sys
import tempfile


def calc_checksum(fname, *, alg_name, chunksize=4096):
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


def download_file(tmpdirname, url, quiet_download):
    if shutil.which("wget") is not None:
        _cmd = ["wget"]
        if quiet_download:
            _cmd.append("--quiet")
    elif shutil.which("curl") is not None:
        _cmd = ["curl", "--remote-name"]
        if quiet_download:
            _cmd.append("--silent")
    else:
        # PR #235 demonstrates how to do this with standard python library
        raise RuntimeError("wget and curl aren't installed")

    _cmd.append(url)
    subprocess.check_output(_cmd, cwd=tmpdirname)
    archive_path = os.path.join(tmpdirname, os.path.basename(url))
    assert os.path.isfile(archive_path)  # sanity check
    return archive_path


def fetch_data(tmpdirname, envvar, url, cksum):
    # retrieve data-dir and make sure it exists
    datadir = os.getenv(envvar, None)
    if datadir is None:
        raise RuntimeError(f"`{envvar}` env var not defined")
    elif not os.path.isabs(datadir):
        raise RuntimeError(f"`{envvar}` env var not an absolute path")
    elif not os.path.isdir(datadir):
        os.mkdir(datadir)

    # download the file
    quiet_download = "true" == os.getenv("CI", None)
    archive_path = download_file(tmpdirname, url, quiet_download)

    # validate checksum
    _measured = calc_checksum(archive_path, alg_name="sha256", chunksize=4096)
    if cksum != _measured:
        raise RuntimeError(
            f"sha1 checksum mismatch\n  actual:   {_measured!r}\n  expected: {cksum!r}"
        )

    subprocess.check_output(["tar", "xzf", archive_path], cwd=datadir)


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as tmpdirname:
        fetch_data(
            tmpdirname,
            "YT_DATA_DIR",
            "https://yt-project.org/data/IsolatedGalaxy.tar.gz",
            "fc081bd4420efd02f7ba2db7eaf4bff0299d5cc4da93436be55a30c219daaf18",
        )
    sys.exit(0)
