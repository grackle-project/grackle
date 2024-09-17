import sys

from .utilities.grdata import main
from .utilities.data_path import _CONFIG_PAIR

if __name__ == '__main__':
    sys.exit(main(*_CONFIG_PAIR, prog_name="python -m pygrackle"))
