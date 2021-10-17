import sys

from ..api import pipeline

import fuc

description = f"""
Run genotyping pipeline for chip data.

Usage examples:
  $ pypgx {fuc.api.common._script_name()} \\
    CYP3A5 \\
    CYP3A5-pipeline \\
    --variants variants.vcf
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Run genotyping pipeline for chip data.',
        description=description,
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
        '--force',
        action='store_true',
        help='Overwrite output directory if it already exists.'
    )

def main(args):
    pipeline.run_chip_pipeline(
        args.gene, args.output, args.variants, force=args.force
    )
