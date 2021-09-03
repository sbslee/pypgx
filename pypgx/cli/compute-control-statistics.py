import sys

from ..api import utils

import fuc
import pysam

description = f"""
This command will compute various statistics for control gene from BAM data.

Input files must be specified with either '--bam' or '--fn'.

Control gene must be specified with either '--gene' or '--region'.

Usage examples:
  $ pypgx {fuc.api.common._script_name()} --fn bam.list --gene VDR
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Compute various statistics for control gene from BAM data.',
        description=description,
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
        '--gene',
        metavar='TEXT',
        help="Control gene (choices: 'EGFR', 'RYR1', 'VDR')."
    )
    parser.add_argument(
        '--region',
        metavar='TEXT',
        help='Custom control region.'
    )
    parser.add_argument(
        '--assembly',
        metavar='TEXT',
        default='GRCh37',
        help="Reference genome assembly (default: 'GRCh37') (choices: 'GRCh37', 'GRCh38')."
    )

def main(args):
    bam_files = []

    if args.bam is None and args.fn is None:
        raise ValueError(
            "Either the 'bam' or 'fn' parameter must be provided.")
    elif args.bam is not None and args.fn is not None:
        raise ValueError(
            "The 'bam' and 'fn' parameters cannot be used together.")
    elif args.bam is not None and args.fn is None:
        if isinstance(args.bam, str):
            bam_files.append(args.bam)
        else:
            bam_files += args.bam
    else:
        bam_files += fuc.api.common.convert_file2list(args.fn)

    df = utils.load_gene_table()

    if args.gene is not None:
        region = df[df.Gene == args.gene][f'{args.assembly}Region'].values[0]

    if all([fuc.api.pybam.has_chr(x) for x in bam_files]):
        region = 'chr' + region

    cf = fuc.api.pycov.CovFrame.from_bam(
        bam=bam_files, region=region, zero=False
    )

    df = cf.df.iloc[:, 2:].describe()

    sys.stdout.write(df.to_csv(sep='\t'))
