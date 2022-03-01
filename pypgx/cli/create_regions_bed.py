import sys

from ..api import utils

import fuc
import pysam

description = f"""
Create a BED file which contains all regions used by PyPGx.
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        help=
"""Create a BED file which contains all regions used by
PyPGx."""
    )
    parser.add_argument(
        '--assembly',
        metavar='TEXT',
        default='GRCh37',
        help=
"""Reference genome assembly (default: 'GRCh37')
(choices: 'GRCh37', 'GRCh38')."""
    )
    parser.add_argument(
        '--add-chr-prefix',
        action='store_true',
        help=
"""Whether to add the 'chr' string in contig names."""
    )
    parser.add_argument(
        '--merge',
        action='store_true',
        help=
"""Whether to merge overlapping intervals (gene names will
be removed too)."""
    )
    parser.add_argument(
        '--sv-genes',
        action='store_true',
        help=
"""Whether to only return genes with SV."""
    )

def main(args):
    bf = utils.create_regions_bed(
        assembly=args.assembly, add_chr_prefix=args.add_chr_prefix,
        merge=args.merge, sv_genes=args.sv_genes
    )
    sys.stdout.write(bf.to_string())
