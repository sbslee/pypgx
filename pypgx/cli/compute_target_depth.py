import tempfile

from ..api import utils

import fuc
import pysam

description = f"""
Compute read depth for the target gene from BAM files.
"""

epilog = f"""
[Example] For the CYP2D6 gene from WGS data:
  $ pypgx {fuc.api.common._script_name()} \\
  CYP2D6 \\
  read-depth.zip \\
  --bam A.bam B.bam

[Example] For the CYP2D6 gene from targeted sequencing data:
  $ pypgx {fuc.api.common._script_name()} \\
  CYP2D6 \\
  read-depth.zip \\
  --fn bam.txt \\
  --bed probes.bed
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        epilog=epilog,
        help='Compute read depth for the target gene from BAM files.',
    )
    parser.add_argument(
        'gene',
        help='Target gene.'
    )
    parser.add_argument(
        'output',
        help='Archive file with the semantic type \n'
             'CovFrame[ReadDepth].'
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
    archive = utils.compute_target_depth(
        args.gene, bam=args.bam, fn=args.fn, assembly=args.assembly,
        bed=args.bed
    )
    archive.to_file(args.output)
