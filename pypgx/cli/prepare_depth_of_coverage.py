import tempfile

from ..api import utils

import fuc
import pysam

description = f"""
###################################################################
# Create TSV file containing read depth for target genes with SV. #
###################################################################

Input BAM files must be specified with either '--bam' or '--fn', but it's an error to use both.

Usage examples:
  $ fuc {fuc.api.common._script_name()} read-depth.tsv --bam A.bam B.bam
  $ fuc {fuc.api.common._script_name()} read-depth.tsv --fn bam.list
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Create TSV file containing read depth for target genes with SV.',
        description=description,
    )
    parser.add_argument(
        'tsv',
        help='TSV file containing read depth.'
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

def main(args):
    cf = utils.prepare_depth_of_coverage(
        bam=args.bam, fn=args.fn, assembly=args.assembly,
    )
    cf.to_file(args.tsv)
