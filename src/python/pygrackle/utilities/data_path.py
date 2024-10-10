########################################################################
#
# Try to get a path to the Grackle input data.
#
#
# Copyright (c) Grackle Development Team. All rights reserved.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

import io
import os
import sys

from pygrackle.__config__ import _is_editable_installation
from pygrackle.grackle_wrapper import get_grackle_version
from pygrackle.utilities.grdata import (
    make_config_objects,
    VersionDataManager,
    _parse_file_registry,
)
from pygrackle.utilities.misc import dirname

# maybe it would be better to export nothing?
__all__ = ["grackle_data_dir"]


def _get_file_registry_contents():
    if _is_editable_installation():
        fname = os.path.join(
            dirname(os.path.abspath(__file__), 2), "file_registry", "file_registry.txt"
        )
        if not os.path.isfile(fname):
            raise RuntimeError(
                "could not find the file_registry.txt in an editable install."
            )
        return fname

    if (sys.version_info.major, sys.version_info.minor) < (3, 9):
        import importlib_resources as resources
    else:
        from importlib import resources
    ref = resources.files("pygrackle.file_registry") / "file_registry.txt"

    contents = ref.read_text(encoding="utf-8")
    return io.StringIO(contents)


def _make_config_pair(grackle_version=None):
    if grackle_version is None:
        grackle_version = get_grackle_version()["version"]
    return make_config_objects(
        grackle_version=grackle_version,
        file_registry_file=_get_file_registry_contents(),
    )


_CONFIG_PAIR = _make_config_pair()
_MANAGER = VersionDataManager.create(*_CONFIG_PAIR)


def _fnames_in_registry():
    # used for testing/debugging
    return tuple(_parse_file_registry(_CONFIG_PAIR[1].file_registry_file).keys())


def _download_all_datafiles():
    """Download all datafiles if it hasn't been downloaded already."""
    registry = _parse_file_registry(_CONFIG_PAIR[1].file_registry_file)
    return _MANAGER.fetch_all(registry)


grackle_data_dir = _MANAGER.version_dir
