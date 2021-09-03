import sys

from ..api import utils

import fuc
import pysam

description = f"""
This command will estimate haplotype phase of observed variants with the Beagle program.

Usage examples:
  $ pypgx {fuc.api.common._script_name()} CYP2D6 in.vcf ref.vcf out
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Estimate haplotype phase of observed variants with the Beagle program.',
        description=description,
    )
    parser.add_argument(
        'gene',
        help='Target gene.'
    )
    parser.add_argument(
        'vcf',
        help='Input VCF file.'
    )
    parser.add_argument(
        'panel',
        help='Reference haplotype panel.'
    )
    parser.add_argument(
        'output',
        help="Output prefix for phased VCF file. For example, '/path/to/output' will generate '/path/to/output.vcf'."
    )
    parser.add_argument(
        '--assembly',
        metavar='TEXT',
        default='GRCh37',
        help="Reference genome assembly (default: 'GRCh37') (choices: 'GRCh37', 'GRCh38')."
    )
    parser.add_argument(
        '--impute',
        action='store_true',
        help='Whether to perform imputation of missing genotypes.'
    )

def main(args):
    utils.estimate_phase_beagle(
        args.gene, args.vcf, args.panel, args.output, assembly=args.assembly,
        impute=args.impute
    )
