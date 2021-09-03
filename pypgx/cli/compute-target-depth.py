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
        'gene',
        help='Target gene.'
    )
    parser.add_argument(
        'output',
        help='Output prefix for zipped TSV file.'
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

    region = utils.get_region(args.gene, assembly=args.assembly)

    cf = fuc.api.pycov.CovFrame.from_bam(
        bam=bam_files, region=f'{prefix}{region}', zero=True
    )

    cf.to_file(f'{args.output}.tsv.gz')
