import sys

from ..api import utils

import fuc
import pysam

description = f"""
Print the main data of specified archive.
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        help=
"""Print the main data of specified archive."""
    )
    parser.add_argument(
        'input',
        help=
"""Input archive file."""
    )

def main(args):
    utils.print_data(args.input)
