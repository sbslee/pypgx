import sys

from ..api import utils

import fuc

description = f"""
Calculate concordance rate between two genotype results.

The command will only use samples that appear in both genotype results.
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Calculate concordance rate between two genotype results.',
        description=description,
    )
    parser.add_argument(
        'first',
        help='First archive file with the semantic type \n'
             'SampleTable[Results].'
    )

    parser.add_argument(
        'second',
        help='Second archive file with the semantic type \n'
             'SampleTable[Results].'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Whether to print the verbose version of output.'
    )

def main(args):
    utils.compare_genotypes(args.first, args.second, verbose=args.verbose)
