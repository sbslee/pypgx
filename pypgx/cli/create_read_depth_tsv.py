import tempfile

from ..api import utils

import fuc
import pysam

description = f"""
#####################################################
# Compute read depth for target gene with BAM data. #
#####################################################

Input BAM files must be specified with either '--bam' or '--fn', but it's an error to use both.

By default, the input data is assumed to be WGS. If it's targeted sequencing, you must provide a BED file with ``bed`` to indicate probed regions.

Usage examples:
  $ fuc {fuc.api.common._script_name()} gene out.zip --bam A.bam B.bam
  $ fuc {fuc.api.common._script_name()} gene out.zip --fn bam.list
  $ fuc {fuc.api.common._script_name()} gene out.zip --fn bam.list --assembly GRCh38
  $ fuc {fuc.api.common._script_name()} gene out.zip --fn bam.list --bed panel.bed
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Compute read depth for target gene with BAM data.',
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
    cf = utils.create_read_depth_tsv(
        bam=args.bam, fn=args.fn, assembly=args.assembly,
    )
    cf.to_file(args.tsv)
