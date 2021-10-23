import sys

from ..api import utils

import fuc

description = f"""
Combine various results for the target gene.
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
        help='Archive file with the semantic type \n'
             'SampleTable[Genotypes].'
    )
    parser.add_argument(
        '--phenotypes',
        metavar='PATH',
        help='Archive file with the semantic type \n'
             'SampleTable[Phenotypes].'
    )
    parser.add_argument(
        '--alleles',
        metavar='PATH',
        help='Archive file with the semantic type \n'
             'SampleTable[Alleles].'
    )
    parser.add_argument(
        '--cnv-calls',
        metavar='PATH',
        help='Archive file with the semantic type \n'
             'SampleTable[CNVCalls].'
    )

def main(args):
    archive = utils.combine_results(
        genotypes=args.genotypes, phenotypes=args.phenotypes,
        alleles=args.alleles, cnv_calls=args.cnv_calls
    )
    archive.to_file(args.results)
