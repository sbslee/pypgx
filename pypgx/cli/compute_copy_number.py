import sys

from ..api import utils

import fuc
import pysam

description = f"""
########################################################
# Compute copy number from read depth for target gene. #
########################################################

The command will convert read depth to copy number by performing intra-sample normalization with control statistics.

If the input data was generated with targeted sequencing, as opposed to WGS, the command will also apply inter-sample normalization using the population statistics. However, for best results it's recommended to manually specify a list of known samples without SV.

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
        help='Archive file with the semantic type CovFrame[ReadDepth].'
    )
    parser.add_argument(
        'control',
        help='Archive file with the semantic type TSV[Statistics].'
    )
    parser.add_argument(
        'output',
        help='Archive file with the semantic type CovFrame[CopyNumber].'
    )
    parser.add_argument(
        '--samples',
        metavar='TEXT',
        nargs='+',
        help='List of known samples with no SV.'
    )

def main(args):
    result = utils.compute_copy_number(
        args.target, args.control, samples=args.samples
    )
    result.to_file(args.output)
