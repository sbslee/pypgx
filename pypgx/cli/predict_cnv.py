import sys

from ..api import utils

import fuc
import pysam

description = f"""
##########################################################
# Predict CNV for target gene based on copy number data. #
##########################################################

If there are missing values because, for example, the input data was generated with targeted sequencing, they will be filled in with the sample's median copy number.

Usage examples:
  $ pypgx {fuc.api.common._script_name()} CYP2D6-copy-number.zip CYP2D6-cnv-calls.zip
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Predict CNV for target gene based on copy number data.',
        description=description,
    )
    parser.add_argument(
        'copy_number',
        metavar='copy-number',
        help='Archive file with the semantic type CovFrame[CopyNumber].'
    )
    parser.add_argument(
        'cnv_calls',
        metavar='cnv-calls',
        help='Archive file with the semantic type SampleTable[CNVCalls].'
    )

def main(args):
    archive = utils.predict_cnv(args.copy_number)
    archive.to_file(args.cnv_calls)
