import sys

from ..api import plot

import fuc
import pysam

description = f"""
###############################################
# Plot allele fraction profile with VCF data. #
###############################################

Usage examples:
  $ pypgx {fuc.api.common._script_name()} CYP2D6 in.vcf
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Plot allele fraction profile with VCF data.',
        description=description,
    )
    parser.add_argument(
        'imported_variants',
        metavar='imported-variants',
        help='VCF file.'
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
    plot.plot_vcf_allele_fraction(
        args.imported_variants, path=args.path, samples=args.samples,
        ymin=args.ymin, ymax=args.ymax
    )
