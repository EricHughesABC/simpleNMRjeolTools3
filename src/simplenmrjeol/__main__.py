"""Enables `python -m simplenmrjeol` as an alternative to the
`simplenmr-jeol` console script installed by pip."""

import sys

from .jason_simpleNMR_cli import main

if __name__ == "__main__":
    sys.exit(main())
