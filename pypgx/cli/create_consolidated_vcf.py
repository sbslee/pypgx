import sys

from ..api import utils

import fuc
import pysam

description = f"""
############################
# Create consolidated VCF. #
############################

Usage examples:
  $ pypgx {fuc.api.common._script_name()} imported.zip phased.zip out.zip
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Create consolidated VCF.',
        description=description,
    )
    parser.add_argument(
        'imported',
        help='Archive file with the semantic type VcfFrame[Imported].'
    )
    parser.add_argument(
        'phased',
        help='Archive file with the semantic type VcfFrame[Phased]'
    )
    parser.add_argument(
        'output',
        help='Archive file with the semantic type VcfFrame[Consolidated].'
    )

def main(args):
    result = utils.create_consolidated_vcf(
        args.imported, args.phased
    )
    result.to_file(args.output)
