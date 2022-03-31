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

def main(args):
    utils.create_input_vcf(
        args.vcf, args.fasta, args.bams, assembly=args.assembly
    )
