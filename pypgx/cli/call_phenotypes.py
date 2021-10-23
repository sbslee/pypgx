import sys

from ..api import utils

import fuc

description = f"""
Call phenotypes for the target gene.
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Call phenotypes for the target gene.',
        description=description,
    )
    parser.add_argument(
        'genotypes',
        help='Archive file with the semantic type SampleTable[Genotypes].'
    )
    parser.add_argument(
        'phenotypes',
        help='Archive file with the semantic type SampleTable[Phenotypes].'
    )

def main(args):
    archive = utils.call_phenotypes(args.genotypes)
    archive.to_file(args.phenotypes)
