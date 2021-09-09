import sys

from ..api import utils

import fuc
import pysam

description = f"""
######################################
# Test a CNV caller for target gene. #
######################################

Usage examples:
  $ pypgx {fuc.api.common._script_name()} caller.zip target.zip calls.zip
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
        help='Archive file with the semantic type SampleTable[CNVCalls].'
    )

def main(args):
    utils.test_cnv_caller(args.caller, args.target, args.calls)
