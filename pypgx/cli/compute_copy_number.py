import sys

from ..api import utils

import fuc
import pysam

description = f"""
Compute copy number from read depth for target gene.

The command will convert read depth to copy number by performing intra-sample
normalization using summary statistics from the control gene.

During copy number analysis, if the input data is targeted sequencing, the
command will apply inter-sample normalization using summary statistics across
all samples. For best results, it is recommended to specify known samples
without SV using --samples-without-sv.
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        help=
"""Compute copy number from read depth for target gene."""
    )
    parser.add_argument(
        'read_depth',
        metavar='read-depth',
        help=
"""Input archive file with the semantic type
CovFrame[ReadDepth]."""
    )
    parser.add_argument(
        'control_statistics',
        metavar='control-statistics',
        help=
"""Input archive file with the semantic type
SampleTable[Statistics]."""
    )
    parser.add_argument(
        'copy_number',
        metavar='copy-number',
        help=
"""Output archive file with the semantic type
CovFrame[CopyNumber]."""
    )
    parser.add_argument(
        '--samples-without-sv',
        metavar='TEXT',
        nargs='+',
        help=
"""List of known samples with no SV."""
    )

def main(args):
    result = utils.compute_copy_number(
        args.read_depth, args.control_statistics,
        samples_without_sv=args.samples_without_sv
    )
    result.to_file(args.copy_number)
