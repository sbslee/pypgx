import sys

from ..api import utils

import fuc

description = f"""
Call SNVs/indels from BAM files for all target genes.

To save computing resources, this method will call variants only for target
genes whose at least one star allele is defined by SNVs/indels. Therefore,
variants will not be called for target genes that have star alleles defined
only by structural variation (e.g. UGT2B17).
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        help=
"""Call SNVs/indels from BAM files for all target genes."""
    )
    parser.add_argument(
        'vcf',
        help=
"""Output VCF file. It must have .vcf.gz as suffix."""
    )
    parser.add_argument(
        'fasta',
        help=
"""Reference FASTA file."""
    )
    parser.add_argument(
        'bams',
        nargs='+',
        help=
"""One or more input BAM files. Alternatively, you can
provide a text file (.txt, .tsv, .csv, or .list)
containing one BAM file per line."""
    )
    parser.add_argument(
        '--assembly',
        metavar='TEXT',
        default='GRCh37',
        help=
"""Reference genome assembly (default: 'GRCh37')
(choices: 'GRCh37', 'GRCh38')."""
    )
    parser.add_argument(
        '--genes',
        metavar='TEXT',
        nargs='+',
        help=
"""List of genes to include."""
    )
    parser.add_argument(
        '--exclude',
        action='store_true',
        help=
"""Exclude specified genes. Ignored when --genes is not
used."""
    )
    parser.add_argument(
        '--dir-path',
        metavar='PATH',
        help=
"""By default, intermediate files (likelihoods.bcf,
calls.bcf, and calls.normalized.bcf) will be stored
in a temporary directory, which is automatically
deleted after creating final VCF. If you provide a
directory path, intermediate files will be stored
there."""
    )
    parser.add_argument(
        '--max-depth',
        metavar='INT',
        type=int,
        default=250,
        help=
"""At a position, read maximally this number of reads
per input file (default: 250). If your input data is
from WGS (e.g. 30X), you don't need to change this
option. However, if it's from targeted sequencing
with ultra-deep coverage (e.g. 500X), then you need
to increase the maximum depth."""
    )

def main(args):
    utils.create_input_vcf(
        args.vcf, args.fasta, args.bams, assembly=args.assembly,
        genes=args.genes, exclude=args.exclude, dir_path=args.dir_path,
        max_depth=args.max_depth
    )
