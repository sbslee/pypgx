import sys

from ..api import utils

import fuc
import pysam

description = f"""
Import read depth data for target gene.
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        help=
"""Import read depth data for target gene."""
    )
    parser.add_argument(
        'gene',
        help=
"""Target gene."""
    )
    parser.add_argument(
        'depth_of_coverage',
        metavar='depth-of-coverage',
        help=
"""Input archive file with the semantic type
CovFrame[DepthOfCoverage]."""
    )
    parser.add_argument(
        'read_depth',
        metavar='read-depth',
        help=
"""Output archive file with the semantic type
CovFrame[ReadDepth]."""
    )
    parser.add_argument(
        '--samples',
        metavar='TEXT',
        nargs='+',
        help=
"""Specify which samples should be included for analysis
by providing a text file (.txt, .tsv, .csv, or .list)
containing one sample per line. Alternatively, you can
provide a list of samples."""
    )
    parser.add_argument(
        '--exclude',
        action='store_true',
        help=
"""Exclude specified samples."""
    )

def main(args):
    archive = utils.import_read_depth(
        args.gene, args.depth_of_coverage, samples=args.samples,
        exclude=args.exclude
    )
    archive.to_file(args.read_depth)
