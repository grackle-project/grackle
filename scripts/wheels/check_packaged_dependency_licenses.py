#!/usr/bin/env python
"""
Check that the vendored dependencies & binary licenses all match

This is intended to be installed once the gracklepy wheel has been installed
"""

from dataclasses import dataclass
import fnmatch
import re
import os
import pathlib
import platform
from typing import Optional
import sys

_IS_LINUX = sys.platform.startswith('linux')
_IS_LINUX_AARCH64 = _IS_LINUX and platform.machine() == "aarch64"
_IS_LINUX_GLIBC = _IS_LINUX and platform.libc_ver()[0] == "glibc"
# note: while platform.libc_ver() can reliably detect whether we're using glibc,
#       it may provide "" rather than the actual name for some other libraries


@dataclass
class LicenseInfo:
    name: str
    location: str  # location where the info was parsed from
    availability: Optional[str] = None
    description: Optional[str] = None
    files: Optional[str] = None # specifies file(s) that the license applies to


def parse_license_file(path):

    def _chunks(f):
        chunk = None
        for lineno, line in enumerate(f, start=1):
            m = re.match(r"^(?P<field>[A-Za-z]+):\s*(?P<value>\S(.*\S)?)$", line)
            if m is None:
                yield chunk
                chunk = None
            elif m.group("field").lower() == "name":
                yield chunk
                chunk = LicenseInfo(name=m.group("value"), location=f"{path!s}:{lineno}")
            elif m is not None:
                setattr(chunk, m.group("field").lower(), m.group("value"))
        yield chunk

    out = {}
    with open(path, "r") as f:
        return [chunk for chunk in _chunks(f) if chunk is not None]


def get_vendored_libs(gracklepy_dir):
    if platform.system() == "Darwin":
        lib_dir = gracklepy_dir / ".dylibs"
    else:
        lib_dir = gracklepy_dir.parent / "gracklepy.libs"
    lib_paths = []
    with os.scandir(lib_dir) as it:
        for entry in it:
            if entry.is_file():
                lib_paths.append(entry.path)
            else:
                raise RuntimeError(f"{lib_dir} holds a non-file: {entry.path}")
    return set(lib_paths)


def main():
    try:
        import gracklepy
    except ImportError:
        print("ERROR: this check requires gracklepy to be installed")
        return 1
    gracklepy_dir = pathlib.Path(gracklepy.__file__).parent

    # get a list of external libraries that we have packaged
    lib_paths = get_vendored_libs(gracklepy_dir)

    # find and parse the license file
    distinfo = list(gracklepy_dir.parent.glob("gracklepy-*.dist-info"))[0]
    license_entries = parse_license_file(distinfo / "licenses" / "LICENSE")

    # now, we will go through and make sure every license_entry corresponds to 1 or more
    # lib_path entries (and vice-versa)

    all_matches = set()
    for entry in license_entries:
        patterns = entry.files.split()  # patterns may be separated by whitespace
        for pat in patterns:
            matches = fnmatch.filter(lib_paths, str(gracklepy_dir.parent / pat))
            if len(matches) > 0:
                all_matches.update(matches)
            elif pat.endswith("libgcc_s*") and _IS_LINUX_GLIBC:
                continue # libgcc_s needs to be distributed for musl, not glibc
            elif pat.endswith("libquadmath*") and _IS_LINUX_AARCH64:
                continue
            else:
                print(
                    f"ERROR: the file pattern `{pat}` doesn't describe a packaged "
                    f"library. This pattern describes `{entry.location}`, which is "
                    f"defined, starting at `{entry.location}`"
                )
                return 1
            all_matches.update(matches)
    if len(all_matches) != len(lib_paths):
        no_matches = lib_paths.difference(all_matches)
        tmp = "ERROR: packaged libs don't have license entries"
        print(tmp, *no_matches, sep="\n   ")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
