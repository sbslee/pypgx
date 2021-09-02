import sys

from ..api import utils

import fuc
import pysam

description = f"""
This command will compute read depth for specified gene from input SAM/BAM/CRAM files.

Usage examples:
  $ fuc {fuc.api.common._script_name()} in.bam
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Compute read depth from SAM/BAM/CRAM files.',
        description=description,
    )
    parser.add_argument(
        'gene',
        help='Target gene.'
    )
    parser.add_argument(
        '--bam',
        metavar='PATH',
        nargs='+',
        help='One or more input files.'
    )
    parser.add_argument(
        '--fn',
        metavar='PATH',
        help='File containing one input filename per line.'
    )
    parser.add_argument(
        '--assembly',
        metavar='TEXT',
        default='GRCh37',
        help="Reference genome assembly (default: 'GRCh37') (choices: 'GRCh37', 'GRCh38')."
    )

def main(args):
    df = utils.load_gene_table()
    s = df[df.Gene == args.gene]
    region = s[f'{args.assembly}Region'].values[0]
    cf = fuc.api.pycov.CovFrame.from_bam(
        bam=args.bam, fn=args.fn, region=region, zero=True
    )
    sys.stdout.write(cf.to_string())
