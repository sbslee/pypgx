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
  $ pypgx {fuc.api.common._script_name()} input.zip output.zip
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Predict CNV for target gene based on copy number data.',
        description=description,
    )
    parser.add_argument(
        'input',
        help='Archive file with the semantic type CovFrame[CopyNumber].'
    )
    parser.add_argument(
        'output',
        help='Archive file with the semantic type SampleTable[CNVCalls].'
    )

def main(args):
    result = utils.predict_cnv(args.input)
    result.to_file(args.output)
