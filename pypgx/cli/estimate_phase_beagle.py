import sys

from ..api import utils

import fuc
import pysam

description = f"""
This command will estimate haplotype phase of observed variants with the Beagle program.

Usage examples:
  $ pypgx {fuc.api.common._script_name()} in.zip ref.vcf out.zip
  $ pypgx {fuc.api.common._script_name()} in.zip ref.vcf out.zip --impute
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Estimate haplotype phase of observed variants with the Beagle program.',
        description=description,
    )
    parser.add_argument(
        'target',
        help='Result file with the semantic type VCF[Imported].'
    )
    parser.add_argument(
        'panel',
        help='Reference haplotype panel.'
    )
    parser.add_argument(
        'output',
        help='Result file with the semantic type VCF[Phased].'
    )
    parser.add_argument(
        '--impute',
        action='store_true',
        help='Whether to perform imputation of missing genotypes.'
    )

def main(args):
    result = utils.estimate_phase_beagle(
        args.target, args.panel, impute=args.impute
    )
    result.to_file(args.output)
