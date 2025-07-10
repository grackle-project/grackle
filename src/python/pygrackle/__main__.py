import sys

from .utilities.grdata import main as grdata_main
from .utilities.data_path import _make_config


def main(args=None):
    return grdata_main(_make_config(), prog_name="python -m pygrackle", args=args)


if __name__ == "__main__":
    sys.exit(main())
