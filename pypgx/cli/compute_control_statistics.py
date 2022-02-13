import sys

from ..api import utils

import fuc
import pysam

description = """
Compute summary statistics for control gene from BAM files.
"""

epilog = f"""
[Example] For the VDR gene from WGS data:
  $ pypgx {fuc.api.common._script_name()} \\
  control-statistcs.zip \\
  1.bam 2.bam \\
  --gene VDR

[Example] For a custom region from targeted sequencing data:
  $ pypgx {fuc.api.common._script_name()} \\
  control-statistcs.zip \\
  bam.list \\
  --region chr1:100-200 \\
  --bed probes.bed
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Compute summary statistics for control gene from BAM
files."""
    )
    parser.add_argument(
        'control_statistics',
        metavar='control-statistics',
        help=
"""Output archive file with the semantic type
SampleTable[Statistics]."""
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
        '--gene',
        metavar='TEXT',
        help=
"""Control gene (recommended choices: 'EGFR', 'RYR1',
'VDR'). Cannot be used with --region."""
    )
    parser.add_argument(
        '--region',
        metavar='TEXT',
        help=
"""Custom region to use as control gene
('chrom:start-end'). Cannot be used with --gene."""
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
to indicate probed regions. Note that the 'chr'
prefix in BED contig names (e.g. 'chr1' vs. '1') will
be automatically added or removed as necessary to
match the BAM contig names."""
    )

def main(args):
    result = utils.compute_control_statistics(
        args.bams, gene=args.gene, region=args.region,
        assembly=args.assembly, bed=args.bed
    )
    result.to_file(args.control_statistics)
