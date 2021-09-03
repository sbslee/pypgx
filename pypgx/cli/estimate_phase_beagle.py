import sys

from ..api import plot

import fuc
import pysam

description = f"""
This command will estimate haplotype phase of observed variants with the Beagle program.

Usage examples:
  $ pypgx {fuc.api.common._script_name()} CYP2D6 target.tsv control.tsv
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
        'depth',
        help='Read depth file for the target gene.'
    )
    parser.add_argument(
        'control',
        help='Summary statistics file for the control gene.'
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
        args.gene, args.depth, args.control, path=args.path,
        samples=args.samples, ymin=args.ymin, ymax=args.ymax
    )
