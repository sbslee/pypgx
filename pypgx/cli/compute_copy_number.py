import sys

from ..api import utils

import fuc
import pysam

description = f"""
########################################################
# Compute copy number from read depth for target gene. #
########################################################

The method will convert read depth from target gene to copy number by performing intra-sample normalization using summary statistics from control gene.

If the input data was generated with targeted sequencing as opposed to WGS, the method will also apply inter-sample normalization using summary statistics across all samples. For best results, it is recommended to manually specify a list of known reference samples that do not have SV.

Usage examples:
  $ pypgx {fuc.api.common._script_name()} CYP2D6-read-depth.zip VDR-control-statistcs.zip CYP2D6-copy-number.zip
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Compute copy number from read depth for target gene.',
        description=description,
    )
    parser.add_argument(
        'read_depth',
        metavar='read-depth',
        help='Archive file with the semantic type CovFrame[ReadDepth].'
    )
    parser.add_argument(
        'control_statistcs',
        metavar='control-statistcs',
        help='Archive file with the semantic type SampleTable[Statistics].'
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
        args.read_depth, args.control_statistcs, samples=args.samples
    )
    result.to_file(args.output)
