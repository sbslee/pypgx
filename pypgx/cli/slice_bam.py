import sys

from ..api import utils

import fuc
import pysam

description = f"""
Slice BAM file for all genes used by PyPGx.
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        help=
"""Slice BAM file for all genes used by PyPGx."""
    )
    parser.add_argument(
        'input',
        help=
"""Input BAM file. It must be already indexed to allow
random access."""
    )
    parser.add_argument(
        'output',
        help=
"""Output BAM file."""
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
        '--genes',
        metavar='TEXT',
        nargs='+',
        help=
"""List of genes to include."""
    )
    parser.add_argument(
        '--exclude',
        action='store_true',
        help=
"""Exclude specified genes. Ignored when --genes is not
used."""
    )

def main(args):
    utils.slice_bam(
        args.input, args.output, assembly=args.assembly, genes=args.genes,
        exclude=args.exclude
    )
