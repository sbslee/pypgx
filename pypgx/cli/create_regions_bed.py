import sys

from ..api import utils

import fuc
import pysam

description = f"""
Create a BED file which contains all regions used by PyPGx.
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        help=
"""Create a BED file which contains all regions used by
PyPGx."""
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
        '--add-chr-prefix',
        action='store_true',
        help=
"""Whether to add the 'chr' string in contig names."""
    )
    parser.add_argument(
        '--merge',
        action='store_true',
        help=
"""Whether to merge overlapping intervals (gene names
will be removed too)."""
    )
    parser.add_argument(
        '--target-genes',
        action='store_true',
        help=
"""Whether to only return target genes, excluding
control genes and paralogs."""
    )
    parser.add_argument(
        '--sv-genes',
        action='store_true',
        help=
"""Whether to only return target genes whose at least
one star allele is defined by structural variation"""
    )
    parser.add_argument(
        '--var-genes',
        action='store_true',
        help=
"""Whether to only return target genes whose at least
one star allele is defined by SNVs/indels."""
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

def main(args):
    bf = utils.create_regions_bed(
        assembly=args.assembly, add_chr_prefix=args.add_chr_prefix,
        merge=args.merge, target_genes=args.target_genes,
        sv_genes=args.sv_genes, var_genes=args.var_genes, genes=args.genes,
        exclude=args.exclude
    )
    sys.stdout.write(bf.to_string())
