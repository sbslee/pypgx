import tempfile

from ..api import utils

import fuc
import pysam

description = f"""
This command will compute read depth for target gene with BAM data.

Input BAM files must be specified with either '--bam' or '--fn', but not both.

Usage examples:
  $ fuc {fuc.api.common._script_name()} gene out.zip --bam A.bam B.bam
  $ fuc {fuc.api.common._script_name()} gene out.zip --fn bam.list
  $ fuc {fuc.api.common._script_name()} gene out.zip --fn bam.list --assembly GRCh38
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Compute read depth for target gene with BAM data.',
        description=description,
    )
    parser.add_argument(
        'gene',
        help='Target gene.'
    )
    parser.add_argument(
        'output',
        help='Archive file with the semantic type CovFrame[ReadDepth].'
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
    result = utils.compute_target_depth(
        args.gene, bam=args.bam, fn=args.fn, assembly=args.assembly
    )
    result.to_file(args.output)
