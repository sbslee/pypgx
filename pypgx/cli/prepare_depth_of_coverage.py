import tempfile

from ..api import utils

import fuc
import pysam

description = """
Prepare a depth of coverage file for all target genes with SV.
"""

epilog = f"""
[Example] When the input data is WGS:
  $ pypgx {fuc.api.common._script_name()} \\
  depth-of-coverage.zip \\
  --bam A.bam B.bam

[Example] When the input data is targeted sequencing:
  $ pypgx {fuc.api.common._script_name()} \\
  depth-of-coverage.zip \\
  --fn bam.txt \\
  --bed probes.bed
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        epilog=epilog,
        help='Prepare a depth of coverage file for all target genes with SV.',
    )
    parser.add_argument(
        'depth_of_coverage',
        metavar='depth-of-coverage',
        help='Archive file with the semantic type \n'
             'CovFrame[DepthOfCoverage].'
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
        help='File containing one BAM file per line. Cannot be used \n'
             'with --bam.'
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
    archive = utils.prepare_depth_of_coverage(
        bam=args.bam, fn=args.fn, assembly=args.assembly, bed=args.bed
    )
    archive.to_file(args.depth_of_coverage)
