import sys

from ..api import genotype

import fuc
import pysam

description = f"""
###################################
# Call genotypes for target gene. #
###################################

Usage examples:
  $ pypgx {fuc.api.common._script_name()} alleles.zip cnv-calls.zip genotypes.zip
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Call genotypes for target gene.',
        description=description,
    )
    parser.add_argument(
        'alleles',
        help='Archive file with the semantic type SampleTable[Alleles].'
    )
    parser.add_argument(
        'cnv_calls',
        help='Archive file with the semantic type SampleTable[CNVCalls].'
    )
    parser.add_argument(
        'genotypes',
        help='Archive file with the semantic type SampleTable[Genotypes].'
    )

def main(args):
    archive = genotype.call_genotypes(alleles=args.alleles, cnv=args.cnv_calls)
    archive.to_file(args.genotypes)
