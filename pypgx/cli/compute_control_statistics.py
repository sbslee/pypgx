import sys

from ..api import utils

import fuc
import pysam

description = f"""
##############################################################
# Compute various statistics for control gene with BAM data. #
##############################################################

Input BAM files must be specified with either '--bam' or '--fn', but it's an error to use both. Similarly, control gene must be specified with either '--gene' or '--region', but it's an error to use both.

By default, the input data is assumed to be WGS. If it's targeted sequencing, you must provide a BED file with '--bed' to indicate probed regions.

Usage examples:
  $ pypgx {fuc.api.common._script_name()} out.zip --bam A.bam B.bam --gene VDR
  $ pypgx {fuc.api.common._script_name()} out.zip --fn bam.list --region chr:start-end
  $ pypgx {fuc.api.common._script_name()} out.zip --fn bam.list --region chr:start-end --assembly GRCh38
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Compute various statistics for control gene with BAM data.',
        description=description,
    )
    parser.add_argument(
        'output',
        help='Archive file with the semantic type SampleTable[Statistics].'
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
        '--gene',
        metavar='TEXT',
        help="Control gene (recommended choices: 'EGFR', 'RYR1', 'VDR')."
    )
    parser.add_argument(
        '--region',
        metavar='TEXT',
        help="Custom region to use as control gene ('chrom:start-end')."
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
    result = utils.compute_control_statistics(
        bam=args.bam, fn=args.fn, gene=args.gene, region=args.region,
        assembly=args.assembly, bed=args.bed
    )
    result.to_file(args.output)
