import sys

from ..api import utils

import fuc
import pysam

description = f"""
Import read depth data for the target gene.
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        help='Import read depth data for the target gene.',
    )
    parser.add_argument(
        'gene',
        help='Target gene.'
    )
    parser.add_argument(
        'depth_of_coverage',
        metavar='depth-of-coverage',
        help='Archive file with the semantic type \n'
             'CovFrame[DepthOfCoverage].'
    )
    parser.add_argument(
        'read_depth',
        metavar='read-depth',
        help='Archive file with the semantic type CovFrame[ReadDepth].'
    )

def main(args):
    archive = utils.import_read_depth(
        args.gene, args.depth_of_coverage
    )
    archive.to_file(args.read_depth)
