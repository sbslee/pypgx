import sys

from ..api import utils

import fuc
import pysam

description = f"""
###########################################
# Import read depth data for target gene. #
###########################################

Usage examples:
  $ pypgx {fuc.api.common._script_name()} CYP2D6 read-depth.tsv CYP2D6-read-depth.zip
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Import read depth data for target gene.',
        description=description,
    )
    parser.add_argument(
        'gene',
        help='Target gene.'
    )
    parser.add_argument(
        'read_depth',
        help='TSV file containing read depth (zipped or unzipped).'
    )
    parser.add_argument(
        'output',
        help='Archive file with the semantic type CovFrame[ReadDepth].'
    )
    parser.add_argument(
        '--assembly',
        metavar='TEXT',
        default='GRCh37',
        choices=['GRCh37', 'GRCh38'],
        help="Reference genome assembly (default: 'GRCh37') (choices: 'GRCh37', 'GRCh38')."
    )
    parser.add_argument(
        '--platform',
        metavar='TEXT',
        default='WGS',
        choices=['WGS', 'Targeted'],
        help="NGS platform (default: 'WGS') (choices: 'WGS', 'Targeted')."
    )

def main(args):
    archive = utils.import_read_depth(
        args.gene, args.read_depth, assembly=args.assembly,
        platform=args.platform
    )
    archive.to_file(args.output)
