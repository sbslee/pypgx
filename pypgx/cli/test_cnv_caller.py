import sys

from ..api import utils

import fuc
import pysam

description = f"""
This command will test a CNV caller for target gene.

Usage examples:
  $ pypgx {fuc.api.common._script_name()} input.zip output.zip
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Test a CNV caller for target gene.',
        description=description,
    )
    parser.add_argument(
        'caller',
        help='Archive file with the semantic type Model[CNV].'
    )
    parser.add_argument(
        'target',
        help='Archive file with the semantic type CovFrame[CopyNumber].'
    )
    parser.add_argument(
        'calls',
        help='Archive file with the semantic type TSV[CNVCalls].'
    )

def main(args):
    utils.test_cnv_caller(args.caller, args.target, args.calls)
