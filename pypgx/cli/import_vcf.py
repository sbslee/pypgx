import sys

from ..api import utils

import fuc
import pysam

description = f"""
####################################
# Import VCF data for target gene. #
####################################

Usage examples:
  $ pypgx {fuc.api.common._script_name()} gene in.vcf out.zip
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Import VCF data for target gene.',
        description=description,
    )
    parser.add_argument(
        'gene',
        help='Target gene.'
    )
    parser.add_argument(
        'vcf',
        help='VCF file (zipped or unzipped).'
    )
    parser.add_argument(
        'output',
        help='Archive file with the semantic type VcfFrame[Imported].'
    )
    parser.add_argument(
        '--assembly',
        metavar='TEXT',
        default='GRCh37',
        help="Reference genome assembly (default: 'GRCh37') (choices: 'GRCh37', 'GRCh38')."
    )

def main(args):
    archive = utils.import_vcf(args.gene, args.vcf, assembly=args.assembly)
    archive.to_file(args.output)    
