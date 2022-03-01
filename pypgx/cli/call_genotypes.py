import sys

from ..api import genotype

import fuc

description = f"""
Call genotypes for target gene.
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        help=
"""Call genotypes for target gene."""
    )
    parser.add_argument(
        'genotypes',
        help=
"""Output archive file with the semantic type
SampleTable[Genotypes]."""
    )
    parser.add_argument(
        '--alleles',
        metavar='PATH',
        help=
"""Input archive file with the semantic type
SampleTable[Alleles]."""
    )
    parser.add_argument(
        '--cnv-calls',
        metavar='PATH',
        help=
"""Input archive file with the semantic type
SampleTable[CNVCalls]."""
    )

def main(args):
    archive = genotype.call_genotypes(
        alleles=args.alleles, cnv_calls=args.cnv_calls
    )
    archive.to_file(args.genotypes)
