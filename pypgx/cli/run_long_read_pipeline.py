import sys

from ..api import pipeline

import fuc

description = f"""
Run genotyping pipeline for long-read sequencing data.
"""

epilog = f"""
[Example] To genotype the CYP3A5 gene from long-read sequencing data:
  $ pypgx {fuc.api.common._script_name()} \\
  CYP3A5 \\
  CYP3A5-pipeline \\
  variants.vcf.gz
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Run genotyping pipeline for long-read sequencing data."""
    )
    parser.add_argument(
        'gene',
        help=
"""Target gene."""
    )
    parser.add_argument(
        'output',
        help=
"""Output directory."""
    )
    parser.add_argument(
        'variants',
        help=
"""Input VCF file must be already BGZF compressed (.gz)
and indexed (.tbi) to allow random access."""
    )
    parser.add_argument(
        '--assembly',
        metavar='TEXT',
        default='GRCh37',
        help=
"""Reference genome assembly (default: 'GRCh37')
(choices: 'GRCh37', 'GRCh38')."""
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help=
"""Overwrite output directory if it already exists."""
    )
    parser.add_argument(
        '--samples',
        metavar='TEXT',
        nargs='+',
        help=
"""Specify which samples should be included for analysis
by providing a text file (.txt, .tsv, .csv, or .list)
containing one sample per line. Alternatively, you
can provide a list of samples."""
    )
    parser.add_argument(
        '--exclude',
        action='store_true',
        help=
"""Exclude specified samples."""
    )

def main(args):
    pipeline.run_long_read_pipeline(
        args.gene, args.output, args.variants, assembly=args.assembly,
        force=args.force, samples=args.samples, exclude=args.exclude
    )
