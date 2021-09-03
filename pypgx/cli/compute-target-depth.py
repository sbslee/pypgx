import sys
import tempfile

from ..api import utils

import fuc
import pysam

description = f"""
This command will compute read depth for target gene with BAM data.

Input files must be specified with either '--bam' or '--fn'.

Usage examples:
  $ fuc {fuc.api.common._script_name()} CYP2D6 --fn bam.list
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Compute read depth for target gene with BAM data.',
        description=description,
    )
    parser.add_argument(
        '--genes',
        metavar='TEXT',
        nargs='+',
        help='Target genes to include.'
    )
    parser.add_argument(
        '--exclude',
        help='Exclude specified genes.'
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

    if all([fuc.api.pybam.has_chr(x) for x in bam_files]):
        prefix = 'chr'
    else:
        prefix = ''

    if args.genes is None:
        genes = utils.list_genes()
    elif args.exclude:
        genes = [x for x in utils.list_genes() if x not in args.genes]
    else:
        genes = args.genes

    with tempfile.TemporaryDirectory() as t:
        with open(f'{t}/temp.bed', 'w') as f:
            for gene in genes:
                region = utils.get_region(gene, assembly=args.assembly)
                chrom, start, end = fuc.common.parse_region(region)
                f.write(f'{prefix}{chrom}\t{start}\t{end}\n')

        cf = fuc.api.pycov.CovFrame.from_bam(
            bam=bam_files, bed=f'{t}/temp.bed', zero=True
        )

    sys.stdout.write(cf.to_string())
