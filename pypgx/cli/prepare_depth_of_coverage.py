import tempfile

from ..api import utils

import fuc
import pysam

description = f"""
##############################################################
# Prepare a depth of coverage file for target genes with SV. #
##############################################################

Input BAM files must be specified with either '--bam' or '--fn', but it's an error to use both.

By default, the input data is assumed to be WGS. If it's targeted sequencing, you must provide a BED file with '--bed' to indicate probed regions.

Usage examples:
  $ pypgx {fuc.api.common._script_name()} depth-of-coverage.tsv --bam A.bam B.bam
  $ pypgx {fuc.api.common._script_name()} depth-of-coverage.tsv --fn bam.list
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Prepare a depth of coverage file for target genes with SV.',
        description=description,
    )
    parser.add_argument(
        'depth_of_coverage',
        metavar='depth-of-coverage',
        help='Depth of coverage file.'
    )
    parser.add_argument(
        '--bam',
        metavar='PATH',
        nargs='+',
        help='One or more BAM files.'
    )
    parser.add_argument(
        '--fn',
        metavar='PATH',
        help='File containing one BAM file per line.'
    )
    parser.add_argument(
        '--assembly',
        metavar='TEXT',
        default='GRCh37',
        help="Reference genome assembly (default: 'GRCh37') (choices: 'GRCh37', 'GRCh38')."
    )
    parser.add_argument(
        '--bed',
        metavar='PATH',
        help='BED file.'
    )

def main(args):
    cf = utils.prepare_depth_of_coverage(
        bam=args.bam, fn=args.fn, assembly=args.assembly, bed=args.bed
    )
    cf.to_file(args.depth_of_coverage)
