import tempfile

from ..api import utils

import fuc
import pysam

description = f"""
This command will compute read depth for target gene with BAM data.

Input files must be specified with either '--bam' or '--fn'.

Usage examples:
  $ fuc {fuc.api.common._script_name()} CYP2D6 --bam A.bam B.bam
  $ fuc {fuc.api.common._script_name()} CYP2D6 --fn bam.list
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
        help='Result file with the semantic type CovFrame[ReadDepth].'
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
