import sys

from ..api import utils

import fuc
import pysam

description = f"""
Predict CNV for the target gene based on copy number data.

Genomic positions that are missing copy number, because for example the input
data is targeted sequencing, will be imputed with forward filling.
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        help='Predict CNV for the target gene based on copy number data.',
    )
    parser.add_argument(
        'copy_number',
        metavar='copy-number',
        help='Archive file with the semantic type CovFrame[CopyNumber].'
    )
    parser.add_argument(
        'cnv_calls',
        metavar='cnv-calls',
        help='Archive file with the semantic type \n'
             'SampleTable[CNVCalls].'
    )
    parser.add_argument(
        '--cnv-caller',
        metavar='PATH',
        help='Archive file with the semantic type Model[CNV]. By \n'
             'default, a pre-trained CNV caller will be used.'
    )

def main(args):
    archive = utils.predict_cnv(args.copy_number, cnv_caller=args.cnv_caller)
    archive.to_file(args.cnv_calls)
