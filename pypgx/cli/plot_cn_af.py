import sys

from ..api import plot

import fuc
import pysam

description = f"""
Plot both copy number profile and allele fraction profile in one figure.
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        help=
"""Plot both copy number profile and allele fraction
profile in one figure."""
    )
    parser.add_argument(
        'copy_number',
        metavar='copy-number',
        help=
"""Input archive file with the semantic type
CovFrame[CopyNumber]."""
    )
    parser.add_argument(
        'imported_variants',
        metavar='imported-variants',
        help=
"""Input archive file with the semantic type
VcfFrame[Imported]."""
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
        default=-0.3,
        help=
"""Y-axis bottom (default: -0.3)."""
    )
    parser.add_argument(
        '--ymax',
        metavar='FLOAT',
        type=float,
        default=6.3,
        help=
"""Y-axis top (default: 6.3)."""
    )
    parser.add_argument(
        '--fontsize',
        metavar='FLOAT',
        type=float,
        default=25,
        help=
"""Text fontsize (default: 25)."""
    )

def main(args):
    plot.plot_cn_af(
        args.copy_number, args.imported_variants, path=args.path,
        samples=args.samples, ymin=args.ymin, ymax=args.ymax,
        fontsize=args.fontsize
    )
