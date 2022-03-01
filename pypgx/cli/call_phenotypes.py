import sys

from ..api import utils

import fuc

description = f"""
Call phenotypes for target gene.
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        help=
"""Call phenotypes for target gene."""
    )
    parser.add_argument(
        'genotypes',
        help=
"""Input archive file with the semantic type
SampleTable[Genotypes]."""
    )
    parser.add_argument(
        'phenotypes',
        help=
"""Output archive file with the semantic type
SampleTable[Phenotypes]."""
    )

def main(args):
    archive = utils.call_phenotypes(args.genotypes)
    archive.to_file(args.phenotypes)
