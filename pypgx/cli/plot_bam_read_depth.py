import sys

from ..api import plot

import fuc
import pysam

description = f"""
Plot read depth profile with BAM data.
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        help=
"""Plot read depth profile with BAM data."""
    )
    parser.add_argument(
        'read_depth',
        metavar='read-depth',
        help=
"""Input archive file with the semantic type
CovFrame[ReadDepth]."""
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
    parser.add_argument(
        '--fontsize',
        metavar='FLOAT',
        type=float,
        default=25,
        help=
"""Text fontsize (default: 25)."""
    )

def main(args):
    plot.plot_bam_read_depth(
        args.read_depth, path=args.path, samples=args.samples,
        ymin=args.ymin, ymax=args.ymax, fontsize=args.fontsize
    )
