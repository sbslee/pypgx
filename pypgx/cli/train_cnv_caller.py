import sys

from ..api import utils

import fuc
import pysam

description = f"""
###########################################
# Train a CNV caller for the target gene. #
###########################################

This command will return a SVM-based multiclass classifier that makes CNV calls using the one-vs-rest stategy.

Usage examples:
  $ pypgx {fuc.api.common._script_name()} CYP2D6-copy-number.zip CYP2D6-cnv-calls.zip CYP2D6-cnv-caller.zip
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Train a CNV caller for the target gene.',
        description=description,
    )
    parser.add_argument(
        'copy_number',
        metavar='copy-number',
        help='Archive file with the semantic type CovFrame[CopyNumber]'
    )
    parser.add_argument(
        'cnv_calls',
        metavar='cnv-calls',
        help='Archive file with the semantic type SampleTable[CNVCalls].'
    )
    parser.add_argument(
        'cnv_caller',
        metavar='cnv-caller',
        help='Archive file with the semantic type Model[CNV].'
    )
    parser.add_argument(
        '--confusion-matrix',
        metavar='PATH',
        help='Write the confusion matrix as a CSV file.'
    )

def main(args):
    result = utils.train_cnv_caller(
        args.copy_number, args.cnv_calls,
        confusion_matrix=args.confusion_matrix
    )
    result.to_file(args.cnv_caller)
