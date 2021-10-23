import sys

from ..api import utils

import fuc
import pysam

description = """
Compute summary statistics for the control gene from BAM files.
"""

epilog = f"""
[Example] To compute summary statistics for the VDR gene from WGS data:
  $ pypgx {fuc.api.common._script_name()} \\
  control-statistcs-VDR.zip \\
  --gene VDR \\
  --bam A.bam B.bam

[Example] For a custom region from targeted sequencing data:
  $ pypgx {fuc.api.common._script_name()} \\
  control-statistcs-VDR.zip \\
  --gene chr1:100-200 \\
  --fn bam.list \\
  --bed probes.bed
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        epilog=epilog,
        help='Compute summary statistics for the control gene from '
             'BAM files.',
    )
    parser.add_argument(
        'control_statistics',
        metavar='control-statistics',
        help='Archive file with the semantic type \n'
             'SampleTable[Statistics].'
    )
    parser.add_argument(
        '--bam',
        metavar='PATH',
        nargs='+',
        help='One or more BAM files. Cannot be used with --fn.'
    )
    parser.add_argument(
        '--fn',
        metavar='PATH',
        help='File containing one BAM file per line. Cannot be \n'
             'used with --bam.'
    )
    parser.add_argument(
        '--gene',
        metavar='TEXT',
        help="Control gene (recommended choices: 'EGFR', 'RYR1', \n"
             "'VDR'). Cannot be used with --region."
    )
    parser.add_argument(
        '--region',
        metavar='TEXT',
        help="Custom region to use as control gene \n"
             "('chrom:start-end'). Cannot be used with --gene."
    )
    parser.add_argument(
        '--assembly',
        metavar='TEXT',
        default='GRCh37',
        help="Reference genome assembly (default: 'GRCh37') \n"
             "(choices: 'GRCh37', 'GRCh38')."
    )
    parser.add_argument(
        '--bed',
        metavar='PATH',
        help="By default, the input data is assumed to be WGS. If it \n"
             "is targeted sequencing, you must provide a BED file to \n"
             "indicate probed regions."
    )

def main(args):
    result = utils.compute_control_statistics(
        bam=args.bam, fn=args.fn, gene=args.gene, region=args.region,
        assembly=args.assembly, bed=args.bed
    )
    result.to_file(args.control_statistics)
