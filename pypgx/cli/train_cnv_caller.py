import sys

from ..api import utils

import fuc
import pysam

description = f"""
#######################################
# Train a CNV caller for target gene. #
#######################################

This command will return a SVM-based multiclass classifier that implements the one-vs-rest stategy.

Usage examples:
  $ pypgx {fuc.api.common._script_name()} target.zip calls.zip output.zip
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Train a CNV caller for target gene.',
        description=description,
    )
    parser.add_argument(
        'target',
        help='Archive file with the semantic type CovFrame[CopyNumber]'
    )
    parser.add_argument(
        'calls',
        help='Archive file with the semantic type SampleTable[CNVCalls].'
    )
    parser.add_argument(
        'output',
        help='Archive file with the semantic type Model[CNV].'
    )

def main(args):
    result = utils.train_cnv_caller(args.target, args.calls)
    result.to_file(args.output)
