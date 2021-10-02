import sys

from ..api import utils

import fuc
import pysam

description = f"""
##########################################################################
# Estimate haplotype phase of observed variants with the Beagle program. #
##########################################################################



Usage examples:
  $ pypgx {fuc.api.common._script_name()} imported-variants.zip phased-variants.zip
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Estimate haplotype phase of observed variants with the Beagle program.',
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
        help='Archive file with the semantic type VcfFrame[Phased].'
    )
    parser.add_argument(
        '--panel',
        metavar='PATH',
        help='Reference haplotype panel. By default, the 1KGP panel is used.'
    )
    parser.add_argument(
        '--impute',
        action='store_true',
        help='Whether to perform imputation of missing genotypes.'
    )

def main(args):
    result = utils.estimate_phase_beagle(
        args.imported_variants, args.panel, impute=args.impute
    )
    result.to_file(args.phased_variants)
