import sys

from ..api import utils

import fuc
import pysam

description = f"""
##############################################
# Filter Archive file for specified samples. #
##############################################

Samples can be specified with either '--samples' or '--fn', but it's an error to use both.

Usage examples:
  $ pypgx {fuc.api.common._script_name()} in.zip out.zip --samples A B C
  $ pypgx {fuc.api.common._script_name()} in.zip out.zip --samples A B C --exclude
  $ pypgx {fuc.api.common._script_name()} in.zip out.zip --fn samples.list
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Filter Archive file for specified samples.',
        description=description,
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
        help='Space-separated list of samples names (the order matters).'
    )
    parser.add_argument(
        '--exclude',
        action='store_true',
        help='Exclude specified samples.'
    )
    parser.add_argument(
        '--fn',
        metavar='PATH',
        help='File containing one sample name per line.'
    )

def main(args):
    archive = utils.filter_samples(
        args.input, samples=args.samples, exclude=args.exclude, fn=args.fn
    )
    archive.to_file(args.output)
