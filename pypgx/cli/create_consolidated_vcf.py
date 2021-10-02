import sys

from ..api import utils

import fuc
import pysam

description = f"""
############################
# Create consolidated VCF. #
############################

Usage examples:
  $ pypgx {fuc.api.common._script_name()} CYP2D6-imported-variants.zip CYP2D6-phased-variants.zip CYP2D6-consolidated-variants.zip
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Create consolidated VCF.',
        description=description,
    )
    parser.add_argument(
        'imported_variants',
        metavar='imported-variants',
        help='Archive file with the semantic type VcfFrame[Imported].'
    )
    parser.add_argument(
        'phased_variants',
        metavar='phased-variants',
        help='Archive file with the semantic type VcfFrame[Phased]'
    )
    parser.add_argument(
        'consolidated_variants',
        metavar='consolidated-variants',
        help='Archive file with the semantic type VcfFrame[Consolidated].'
    )

def main(args):
    archive = utils.create_consolidated_vcf(
        args.imported_variants, args.phased_variants
    )
    archive.to_file(args.consolidated_variants)
