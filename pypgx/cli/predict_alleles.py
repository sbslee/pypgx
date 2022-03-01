import sys

from ..api import utils

import fuc
import pysam

description = f"""
Predict candidate star alleles based on observed variants.
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        help=
"""Predict candidate star alleles based on observed
variants."""
    )
    parser.add_argument(
        'consolidated_variants',
        metavar='consolidated-variants',
        help=
"""Input archive file with the semantic type
VcfFrame[Consolidated]."""
    )
    parser.add_argument(
        'alleles',
        help=
"""Output archive file with the semantic type
SampleTable[Alleles]."""
    )

def main(args):
    alleles = utils.predict_alleles(args.consolidated_variants)
    alleles.to_file(args.alleles)
