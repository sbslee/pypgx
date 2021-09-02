import sys

from ..api import utils

import fuc
import pysam

description = f"""
This command will compute various statistics for control gene from input SAM/BAM/CRAM files.

Usage examples:
  $ fuc {fuc.api.common._script_name()} in.bam
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Compute various statistics for control gene from SAM/BAM/CRAM files.',
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

def main(args):
    utils.plot_copy_number(
        args.gene, args.depth, args.control, path=args.path, samples=args.samples
    )
