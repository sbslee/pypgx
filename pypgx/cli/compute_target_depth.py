import tempfile

from ..api import utils

import fuc
import pysam

description = f"""
Compute read depth for target gene from BAM files.
"""

epilog = f"""
[Example] For the CYP2D6 gene from WGS data:
  $ pypgx {fuc.api.common._script_name()} \\
  CYP2D6 \\
  read-depth.zip \\
  1.bam 2.bam

[Example] For the CYP2D6 gene from targeted sequencing data:
  $ pypgx {fuc.api.common._script_name()} \\
  CYP2D6 \\
  read-depth.zip \\
  bam.list \\
  --bed probes.bed
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Compute read depth for target gene from BAM files."""
    )
    parser.add_argument(
        'gene',
        help=
"""Target gene."""
    )
    parser.add_argument(
        'read_depth',
        metavar='read-depth',
        help=
"""Output archive file with the semantic type
CovFrame[ReadDepth]."""
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
        '--bed',
        metavar='PATH',
        help=
"""By default, the input data is assumed to be WGS. If it
is targeted sequencing, you must provide a BED file to
indicate probed regions."""
    )

def main(args):
    archive = utils.compute_target_depth(
        args.gene, args.bams, assembly=args.assembly, bed=args.bed
    )
    archive.to_file(args.read_depth)
