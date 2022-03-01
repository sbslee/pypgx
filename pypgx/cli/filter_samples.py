import sys

from ..api import utils

import fuc
import pysam

description = f"""
Filter Archive file for specified samples.
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        help=
"""Filter Archive file for specified samples."""
    )
    parser.add_argument(
        'input',
        help=
"""Input archive file."""
    )
    parser.add_argument(
        'output',
        help=
"""Output archive file."""
    )
    parser.add_argument(
        'samples',
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
    archive = utils.filter_samples(
        args.input, args.samples, exclude=args.exclude
    )
    archive.to_file(args.output)
