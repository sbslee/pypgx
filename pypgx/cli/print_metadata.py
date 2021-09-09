import sys

from ..api import utils

import fuc
import pysam

description = f"""
############################################
# Print the metadata of specified archive. #
############################################

Usage examples:
  $ pypgx {fuc.api.common._script_name()} in.zip
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Print the metadata of specified archive.',
        description=description,
    )
    parser.add_argument(
        'input',
        help='Archive file.'
    )

def main(args):
    utils.print_metadata(args.input)
