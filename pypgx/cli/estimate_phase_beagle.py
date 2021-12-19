import sys

from ..api import utils

import fuc
import pysam

description = f"""
Estimate haplotype phase of observed variants with the Beagle program.
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        help='Estimate haplotype phase of observed variants with the '
             'Beagle program.',
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
        help='VCF file corresponding to a reference haplotype panel \n'
             '(compressed or uncompressed). By default, the 1KGP panel \n'
             'is used.'
    )
    parser.add_argument(
        '--impute',
        action='store_true',
        help='Perform imputation of missing genotypes.'
    )

def main(args):
    result = utils.estimate_phase_beagle(
        args.imported_variants, args.panel, impute=args.impute
    )
    result.to_file(args.phased_variants)
