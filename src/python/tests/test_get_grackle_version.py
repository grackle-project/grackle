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

import os
import subprocess

def query_grackle_version_props():
    # retrieve the current version information with git
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # get the name of the most recent tag preceeding this commit:
    _rslt = subprocess.run(["git", "describe", "--abbrev=0", "--tags"],
                           check = True, capture_output = True)
    most_recent_tag = _rslt.stdout.decode().rstrip()
    if 'grackle-' == most_recent_tag[:8]:
        latest_tagged_version = most_recent_tag[8:]
    else:
        raise RuntimeError(
            "expected the first 8 characters of the most recent git-tag to be "
            "equal to 'grackle-'"
        )

    # get the actual revision when most_recent tag was introduced
    _rslt = subprocess.run(["git", "rev-parse", "-n", "1", most_recent_tag],
                           check = True, capture_output = True)
    revision_of_tag = _rslt.stdout.decode().rstrip()

    # get the branch name and current revision
    _rslt = subprocess.run(["git", "rev-parse", "--abbrev-ref", "HEAD"],
                           check = True, capture_output = True)
    branch = _rslt.stdout.decode().rstrip()

    _rslt = subprocess.run(["git", "rev-parse", "HEAD"],
                           check = True, capture_output = True)
    revision = _rslt.stdout.decode().rstrip()

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
        raise ValueError(
            "get_grackle_version should return a dictionary with the 3 items"
        )
    elif results['branch'] != branch:
        raise ValueError(
            f"expected get_grackle_version()['branch'] to be '{branch}', not "
            f"'{results['branch']}'"
        )
    elif results['revision'] != revision:
        raise ValueError(
            f"expected get_grackle_version()['revision'] to be '{revision}', "
            f"not '{results['revision']}'"
        )

    # now compare version numbers (currently we just check that both versions
    # are identical if the git-tag was updated in the most recent release. in
    # all other cases we simply require that the versions are different)
    if tagged_on_current_revision:
        if results['version'] != latest_tagged_version:
            raise ValueError(
                "expected get_grackle_version()['version'] to be "
                "'{latest_tagged_version}', not '{results['version']}'"
            )
    elif results['version'] == latest_tagged_version:
        raise ValueError(
            f"expected get_grackle_version()['version'] to be different from "
            f"'{latest_tagged_version}', since that version is associated with "
            "a different commit."
        )
