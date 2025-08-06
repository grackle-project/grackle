"""
This is a compatibility layer used to help people to transition to importing
from gracklepy
"""

import warnings
import sys

warnings.warn(
    "pygrackle has been renamed to gracklepy. Please import directly from "
    "gracklepy. In a future release, you will no longer be able to import from "
    "pygrackle.",
    category=DeprecationWarning,
    stacklevel=2,
)

from gracklepy import *  # noqa: E402, F401, F403

# this is hack I read about online...
sys.modules[__name__] = sys.modules["gracklepy"]
