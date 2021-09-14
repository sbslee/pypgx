import sys

from ..api import utils

import fuc
import pysam

description = f"""
##############################################################
# Predict candidate star alleles based on observed variants. #
##############################################################

Usage examples:
  $ pypgx {fuc.api.common._script_name()} CYP2D6-consolidated-variants.zip CYP2D6-alleles.zip
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Predict candidate star alleles based on observed variants.',
        description=description,
    )
    parser.add_argument(
        'consolidated_variants',
        metavar='consolidated-variants',
        help='Archive file with the semantic type VcfFrame[Consolidated].'
    )
    parser.add_argument(
        'alleles',
        help='Archive file with the semantic type SampleTable[Alleles].'
    )

def main(args):
    alleles = utils.predict_alleles(args.consolidated_variants)
    alleles.to_file(args.alleles)
