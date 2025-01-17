import sys

from .utilities.grdata import main as grdata_main
from .utilities.data_path import _make_config_pair

def main(args=None):
    return grdata_main(
        *_make_config_pair(), prog_name="python -m pygrackle", args=args
    )


if __name__ == '__main__':
    sys.exit(main())
