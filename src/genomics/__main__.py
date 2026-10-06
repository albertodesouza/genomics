"""Run the common genomics CLI with `python -m genomics`."""

import sys

from .cli import cli_main


if __name__ == "__main__":
    sys.exit(cli_main())
