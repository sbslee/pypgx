import sys

from ..api import utils

import fuc
import pysam

description = f"""
This command will compute copy number from read depth for target gene.

Usage examples:
  $ pypgx {fuc.api.common._script_name()} target.zip control.zip output.zip
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Compute copy number from read depth for target gene.',
        description=description,
    )
    parser.add_argument(
        'target',
        help='Result file with the semantic type CovFrame[ReadDepth].'
    )
    parser.add_argument(
        'control',
        help='Result file with the semantic type TSV[Statistics].'
    )
    parser.add_argument(
        'output',
        help='Result file with the semantic type CovFrame[CopyNumber].'
    )

def main(args):
    result = utils.compute_copy_number(args.target, args.control)
    result.to_file(args.output)
