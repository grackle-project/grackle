########################################################################
#
# Volumetric heating rate tests
#
#
# Copyright (c) Grackle Development Team. All rights reserved.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

from pygrackle.grackle_wrapper import get_grackle_version
from packaging.version import Version, InvalidVersion

import os
import subprocess

def query_grackle_version_props():
    # retrieve the current version information with git

    cur_dir = os.path.dirname(os.path.abspath(__file__))

    def _run(args):
        _rslt = subprocess.run(
            args, cwd=cur_dir, check=True, capture_output=True,
        )
        return _rslt.stdout.decode().rstrip()

    # get the name of the most recent tag preceeding this commit:
    most_recent_tag = _run(["git", "describe", "--abbrev=0", "--tags"])
    if most_recent_tag.startswith("grackle-"):
        latest_tagged_version = most_recent_tag[8:]
    elif most_recent_tag.startswith("gold-standard-v"):
        latest_tagged_version = most_recent_tag[15:]
    else:
        raise RuntimeError(
            "expected the most recent git-tag to start with "
            "'grackle-' or 'gold-standard-v'"
        )

    # get the actual revision when most_recent tag was introduced
    revision_of_tag = _run(["git", "rev-parse", "-n", "1", most_recent_tag])

    # get the branch name and current revision
    branch = _run(["git", "rev-parse", "--abbrev-ref", "HEAD"])

    revision = _run(["git", "rev-parse", "HEAD"])

    tagged_on_current_revision = revision == revision_of_tag
    return latest_tagged_version, branch, revision, tagged_on_current_revision

def test_get_grackle_version():
    # this test assumes that Grackle was compiled with the currently checked
    # out version of the repository

    latest_tagged_version, branch, revision, tagged_on_current_revision \
        = query_grackle_version_props()

    # perform the easy checks
    results = get_grackle_version()
    if (len(results) != 3):
        raise RuntimeError(
            "get_grackle_version should return a dictionary with the 3 items"
        )
    elif results["branch"] == "N/A" and results["revision"] == "N/A":
        # in this scenario the core library was compiled outside of a git repository
        # so we skip checks of branch and revision
        pass
    elif results['branch'] != branch:
        raise RuntimeError(
            f"expected get_grackle_version()['branch'] to be '{branch}', not "
            f"'{results['branch']}'. Was a different branch of the repository "
            "checked out when the c library was compiled?"
        )
    elif results['revision'] != revision:
        raise RuntimeError(
            f"expected get_grackle_version()['revision'] to be '{revision}', "
            f"not '{results['revision']}'. Was a different revision of the "
            "repository checked out when the c library was compiled?"
        )

    # Check the formatting of the version strings (given by the tag of the most
    # recently tagged commit and given by get_grackle_version()['version'])

    description_string_pairs = [
        ("get_grackle_version()['version']", results['version']),
        ("the tag of the most recently tagged commit", latest_tagged_version)
    ]

    version_objects = []
    for description, str_val in description_string_pairs:
        try:
            # parse is fairly permissive about the file format
            version_objects.append(Version(str_val))
        except InvalidVersion:
            raise RuntimeError(
                f"The value given by {description}, {str_val!r}, doesn't have "
                "the expected formatting."
            )
        # Could add checks enforce other conventions here

    rslt_version_obj, latest_tagged_version_obj = version_objects

    # now check the consistency between the version numbers
    if tagged_on_current_revision:
        if results['version'] != latest_tagged_version:
            raise RuntimeError(
                "expected get_grackle_version()['version'] to be "
                f"{latest_tagged_version!r}, not {results['version']!r}.\n\n"
                "Did you tag a commit with a new version number & forget to "
                "update the version used by the Makefile?"
            )
    elif latest_tagged_version_obj >= rslt_version_obj:
        raise RuntimeError(
            "Something is wrong. The version number given by "
            f"get_grackle_version()['version'], {results['version']!r}, must "
            "exceed the version number encoded in the tag of the most recently "
            f"tagged git commit, {latest_tagged_version!r}.\n\n"
            "Did you update the version used by the Makefile & forget to "
            "update the commit's tag?"
        )
