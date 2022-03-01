import sys

from ..api import plot

import fuc
import pysam

description = f"""
Plot copy number profile from CovFrame[CopyNumber].
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        help=
"""Plot copy number profile from CovFrame[CopyNumber]."""
    )
    parser.add_argument(
        'copy_number',
        metavar='copy-number',
        help=
"""Input archive file with the semantic type
CovFrame[CopyNumber]."""
    )
    parser.add_argument(
        '--fitted',
        action='store_true',
        help=
"""Show the fitted line as well."""
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
    plot.plot_bam_copy_number(
        args.copy_number, fitted=args.fitted, path=args.path,
        samples=args.samples, ymin=args.ymin, ymax=args.ymax,
        fontsize=args.fontsize
    )
