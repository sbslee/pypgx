import sys

from ..api import utils

import fuc
import pysam

description = f"""
Import variant (SNV/indel) data for the target gene.

The command will first slice input VCF for the target gene and then assess
whether every genotype call in the sliced VCF is haplotype phased. It will
return an archive file with the semantic type VcfFrame[Consolidated] if the
VCF is fully phased or otherwise VcfFrame[Imported].
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        help='Import variant (SNV/indel) data for the target gene',
    )
    parser.add_argument(
        'gene',
        help='Target gene.'
    )
    parser.add_argument(
        'vcf',
        help='Input VCF file must be already BGZF compressed (.gz) and \n'
             'indexed (.tbi) to allow random access.'
    )
    parser.add_argument(
        'imported_variants',
        metavar='imported-variants',
        help='Archive file with the semantic type VcfFrame[Imported] \n'
             'or VcfFrame[Consolidated].'
    )
    parser.add_argument(
        '--assembly',
        metavar='TEXT',
        default='GRCh37',
        help="Reference genome assembly (default: 'GRCh37') (choices: \n"
             "'GRCh37', 'GRCh38')."
    )
    parser.add_argument(
        '--platform',
        metavar='TEXT',
        default='WGS',
        choices=['WGS', 'Targeted', 'Chip'],
        help="Genotyping platform (default: 'WGS') (choices: 'WGS', \n"
             "'Targeted', 'Chip')."
    )
    parser.add_argument(
        '--samples',
        metavar='PATH',
        nargs='+',
        help='Specify which samples should be included for analysis \n'
             'by providing a text file (.txt, .tsv, .csv, or .list) \n'
             'containing one sample per line. Alternatively, you can \n'
             'provide a list of samples.'
    )
    parser.add_argument(
        '--exclude',
        action='store_true',
        help='Exclude specified samples.'
    )

def main(args):
    archive = utils.import_variants(
        args.gene, args.vcf, assembly=args.assembly, platform=args.platform,
        samples=args.samples, exclude=args.exclude
    )
    archive.to_file(args.imported_variants)
