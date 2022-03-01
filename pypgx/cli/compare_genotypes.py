import sys

from ..api import utils

import fuc

description = f"""
Calculate concordance between two genotype results.

Only samples that appear in both genotype results will be used to calculate
concordance for genotype calls as well as CNV calls.
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        help=
"""Calculate concordance between two genotype results."""
    )
    parser.add_argument(
        'first',
        help=
"""First archive file with the semantic type
SampleTable[Results]."""
    )

    parser.add_argument(
        'second',
        help=
"""Second archive file with the semantic type
SampleTable[Results]."""
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help=
"""Whether to print the verbose version of output, including
discordant calls."""
    )

def main(args):
    utils.compare_genotypes(args.first, args.second, verbose=args.verbose)
