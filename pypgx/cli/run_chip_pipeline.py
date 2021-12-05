import sys

from ..api import pipeline

import fuc

description = f"""
Run PyPGx's genotyping pipeline for chip data.
"""

epilog = f"""
[Example] To genotype the CYP3A5 gene from chip data:
  $ pypgx {fuc.api.common._script_name()} \\
  CYP3A5 \\
  CYP3A5-pipeline \\
  variants.vcf
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        epilog=epilog,
        help="Run PyPGx's genotyping pipeline for chip data.",
    )
    parser.add_argument(
        'gene',
        help='Target gene.'
    )
    parser.add_argument(
        'output',
        help='Output directory.'
    )
    parser.add_argument(
        'variants',
        help='VCF file (zipped or unzipped).'
    )
    parser.add_argument(
        '--assembly',
        metavar='TEXT',
        default='GRCh37',
        help="Reference genome assembly (default: 'GRCh37') (choices: \n"
             "'GRCh37', 'GRCh38')."
    )
    parser.add_argument(
        '--impute',
        action='store_true',
        help='Perform imputation of missing genotypes.'
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help='Overwrite output directory if it already exists.'
    )

def main(args):
    pipeline.run_chip_pipeline(
        args.gene, args.output, args.variants, assembly=args.assembly,
        impute=args.impute, force=args.force
    )
