import tempfile

from ..api import utils

import fuc
import pysam

description = """
Prepare a depth of coverage file for all target genes with SV from BAM files.
"""

epilog = f"""
[Example] From WGS data:
  $ pypgx {fuc.api.common._script_name()} \\
  depth-of-coverage.zip \\
  1.bam 2.bam

[Example] From targeted sequencing data:
  $ pypgx {fuc.api.common._script_name()} \\
  depth-of-coverage.zip \\
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
"""Prepare a depth of coverage file for all target
genes with SV from BAM files."""
    )
    parser.add_argument(
        'depth_of_coverage',
        metavar='depth-of-coverage',
        help=
"""Output archive file with the semantic type
CovFrame[DepthOfCoverage]."""
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
"""By default, the input data is assumed to be WGS. If
it's targeted sequencing, you must provide a BED file
to indicate probed regions. Note that the 'chr' prefix
in contig names (e.g. 'chr1' vs. '1') will be
automatically added or removed as necessary to match
the input BAM's contig names."""
    )

def main(args):
    archive = utils.prepare_depth_of_coverage(
        args.bams, assembly=args.assembly, bed=args.bed
    )
    archive.to_file(args.depth_of_coverage)
