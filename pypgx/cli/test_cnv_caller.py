import sys

from ..api import utils

import fuc
import pysam

description = f"""
Test a CNV caller for the target gene.
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        help='Test a CNV caller for the target gene.',
    )
    parser.add_argument(
        'cnv_caller',
        metavar='cnv-caller',
        help='Archive file with the semantic type Model[CNV].'
    )
    parser.add_argument(
        'copy_number',
        metavar='copy-number',
        help='Archive file with the semantic type \n'
             'CovFrame[CopyNumber].'
    )
    parser.add_argument(
        'cnv_calls',
        metavar='cnv-calls',
        help='Archive file with the semantic type \n'
             'SampleTable[CNVCalls].'
    )
    parser.add_argument(
        '--confusion-matrix',
        metavar='PATH',
        help='Write the confusion matrix as a CSV file.'
    )

def main(args):
    utils.test_cnv_caller(
        args.cnv_caller, args.copy_number, args.cnv_calls,
        confusion_matrix=args.confusion_matrix
    )
