import sys

from ..api import plot

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
        'vcf',
        help='VCF file.'
    )
    parser.add_argument(
        '--assembly',
        metavar='TEXT',
        default='GRCh37',
        help="Reference genome assembly (default: 'GRCh37') (choices: 'GRCh37', 'GRCh38')."
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
    plot.plot_read_depth(
        args.gene, args.vcf, assembly=args.assembly, path=args.path,
        samples=args.samples, ymin=args.ymin, ymax=args.ymax
    )
