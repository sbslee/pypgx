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
        help='Create a consolidated VCF file.',
        description=description,
    )
    parser.add_argument(
        'imported_variants',
        metavar='imported-variants',
        help='Archive file with the semantic type \n'
             'VcfFrame[Imported].'
    )
    parser.add_argument(
        'phased_variants',
        metavar='phased-variants',
        help='Archive file with the semantic type \n'
             'VcfFrame[Phased].'
    )
    parser.add_argument(
        'consolidated_variants',
        metavar='consolidated-variants',
        help='Archive file with the semantic type \n'
             'VcfFrame[Consolidated].'
    )

def main(args):
    archive = utils.create_consolidated_vcf(
        args.imported_variants, args.phased_variants
    )
    archive.to_file(args.consolidated_variants)
