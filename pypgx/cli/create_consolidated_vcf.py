import sys

from ..api import utils

import fuc
import pysam

description = f"""
Create a consolidated VCF file.
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        help=
"""Create a consolidated VCF file."""
    )
    parser.add_argument(
        'imported_variants',
        metavar='imported-variants',
        help=
"""Input archive file with the semantic type
VcfFrame[Imported]."""
    )
    parser.add_argument(
        'phased_variants',
        metavar='phased-variants',
        help=
"""Input archive file with the semantic type
VcfFrame[Phased]."""
    )
    parser.add_argument(
        'consolidated_variants',
        metavar='consolidated-variants',
        help=
"""Output archive file with the semantic type
VcfFrame[Consolidated]."""
    )

def main(args):
    archive = utils.create_consolidated_vcf(
        args.imported_variants, args.phased_variants
    )
    archive.to_file(args.consolidated_variants)
