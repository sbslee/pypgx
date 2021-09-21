import sys

from ..api import plot

import fuc
import pysam

description = f"""
###########################################
# Plot copy number profile with BAM data. #
###########################################

Usage examples:
  $ pypgx {fuc.api.common._script_name()} CYP2D6-copy-number.zip
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Plot copy number profile with BAM data.',
        description=description,
    )
    parser.add_argument(
        'copy_number',
        metavar='copy-number',
        help='Archive file with the semantic type CovFrame[CopyNumber].'
    )
    parser.add_argument(
        '--path',
        metavar='PATH',
        help='Create plots in this directory.'
    )
    parser.add_argument(
        '--samples',
        metavar='TEXT',
        nargs='+',
        help='Create plots only for these samples.'
    )
    parser.add_argument(
        '--ymin',
        metavar='FLOAT',
        type=float,
        help='Y-axis bottom.'
    )
    parser.add_argument(
        '--ymax',
        metavar='FLOAT',
        type=float,
        help='Y-axis top.'
    )

def main(args):
    plot.plot_bam_copy_number(
        args.copy_number, path=args.path, samples=args.samples,
        ymin=args.ymin, ymax=args.ymax
    )
