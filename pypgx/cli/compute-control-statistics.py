import sys

from ..api import utils

import fuc
import pysam

description = f"""
This command will compute various statistics for control gene with BAM data.

Input files must be specified with either '--bam' or '--fn'.

Control gene must be specified with either '--gene' or '--region'.

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
        help='Archive file with the semantic type TSV[Statistics].'
    )
    parser.add_argument(
        '--bam',
        metavar='PATH',
        nargs='+',
        help='One or more input files.'
    )
    parser.add_argument(
        '--fn',
        metavar='PATH',
        help='File containing one input filename per line.'
    )
    parser.add_argument(
        '--gene',
        metavar='TEXT',
        help="Control gene (choices: 'EGFR', 'RYR1', 'VDR')."
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

def main(args):
    result = utils.compute_control_statistics(
        bam=args.bam, fn=args.fn, gene=args.gene, region=args.region,
        assembly=args.assembly
    )
    result.to_file(args.output)
