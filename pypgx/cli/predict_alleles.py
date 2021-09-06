import sys

from ..api import utils

import fuc
import pysam

description = f"""
This command will predict candidate star alleles based on observed variants.

Usage examples:
  $ pypgx {fuc.api.common._script_name()} input.zip output.zip
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Predict candidate star alleles based on observed variants.',
        description=description,
    )
    parser.add_argument(
        'input',
        help='Result file with the semantic type VcfFrame[Consolidated].'
    )
    parser.add_argument(
        'output',
        help='Result file with the semantic type TSV[Alleles].'
    )

def main(args):
    result = utils.predict_alleles(args.input)
    result.to_file(args.output)
