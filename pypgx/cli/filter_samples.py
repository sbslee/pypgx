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
        help='Filter Archive file for specified samples.',
    )
    parser.add_argument(
        'input',
        help='Input archive file.'
    )
    parser.add_argument(
        'output',
        help='Output archive file.'
    )
    parser.add_argument(
        '--samples',
        metavar='TEXT',
        nargs='+',
        help='List of samples names (the order matters). Cannot be \n'
             'used with --fn.'
    )
    parser.add_argument(
        '--fn',
        metavar='PATH',
        help='File containing one sample name per line. Cannot be \n'
             'used with --samples.'
    )
    parser.add_argument(
        '--exclude',
        action='store_true',
        help='Exclude specified samples.'
    )

def main(args):
    archive = utils.filter_samples(
        args.input, samples=args.samples, fn=args.fn, exclude=args.exclude
    )
    archive.to_file(args.output)
