import sys

from ..api import plot

import fuc
import pysam

description = f"""
Plot read depth profile with VCF data.
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        help=
"""Plot read depth profile with VCF data."""
    )
    parser.add_argument(
        'gene',
        help=
"""Target gene."""
    )
    parser.add_argument(
        'vcf',
        help=
"""Input VCF file."""
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
        '--path',
        metavar='PATH',
        help=
"""Create plots in this directory."""
    )
    parser.add_argument(
        '--samples',
        metavar='TEXT',
        nargs='+',
        help=
"""Specify which samples should be included for analysis
by providing a text file (.txt, .tsv, .csv, or .list)
containing one sample per line. Alternatively, you can
provide a list of samples."""
    )
    parser.add_argument(
        '--ymin',
        metavar='FLOAT',
        type=float,
        help=
"""Y-axis bottom."""
    )
    parser.add_argument(
        '--ymax',
        metavar='FLOAT',
        type=float,
        help=
"""Y-axis top."""
    )

def main(args):
    plot.plot_vcf_read_depth(
        args.gene, args.vcf, assembly=args.assembly, path=args.path,
        samples=args.samples, ymin=args.ymin, ymax=args.ymax
    )
