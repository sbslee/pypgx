import sys

from ..api import utils

import fuc
import pysam

description = f"""
This command will create a BED file which contains all regions used by PyPGx.

Usage examples:
  $ pypgx {fuc.api.common._script_name()} > regions.bed
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Create a BED file which contains all regions used by PyPGx.',
        description=description,
    )
    parser.add_argument(
        '--assembly',
        metavar='TEXT',
        default='GRCh37',
        help="Reference genome assembly (default: 'GRCh37') (choices: 'GRCh37', 'GRCh38')."
    )
    parser.add_argument(
        '--chr_prefix',
        action='store_true',
        help="Whether to add the 'chr' string in contig names."
    )

def main(args):
    bf = utils.create_regions_bed(
        assembly=args.assembly, chr_prefix=args.chr_prefix
    )
    sys.stdout.write(bf.to_string())
