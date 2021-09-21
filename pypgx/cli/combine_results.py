import sys

from ..api import utils

import fuc

description = f"""
################################################
# Combine various results for the target gene. #
################################################

Usage examples:
  $ pypgx {fuc.api.common._script_name()} CYP2D6-results.zip --genotypes CYP2D6-genotypes.zip --alleles CYP2D6-alleles.zip --cnv-calls CYP2D6-cnv-calls.zip
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Combine various results for the target gene.',
        description=description,
    )
    parser.add_argument(
        'results',
        help='Archive file with the semantic type SampleTable[Results].'
    )
    parser.add_argument(
        '--genotypes',
        metavar='PATH',
        help='Archive file with the semantic type SampleTable[Genotypes].'
    )
    parser.add_argument(
        '--alleles',
        metavar='PATH',
        help='Archive file with the semantic type SampleTable[Alleles].'
    )
    parser.add_argument(
        '--cnv-calls',
        metavar='PATH',
        help='Archive file with the semantic type SampleTable[CNVCalls].'
    )

def main(args):
    archive = utils.combine_results(
        genotypes=args.genotypes, alleles=args.alleles,
        cnv_calls=args.cnv_calls
    )
    archive.to_file(args.results)
