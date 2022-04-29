import sys

from ..api import utils

import fuc
import pysam

description = """
Compute summary statistics for control gene from BAM files.

Note that for the arguments gene and --bed, the 'chr' prefix in contig names
(e.g. 'chr1' vs. '1') will be automatically added or removed as necessary to
match the input BAM's contig names.
"""

epilog = f"""
[Example] For the VDR gene from WGS data:
  $ pypgx {fuc.api.common._script_name()} \\
  VDR \\
  control-statistics.zip \\
  1.bam 2.bam

[Example] For a custom region from targeted sequencing data:
  $ pypgx {fuc.api.common._script_name()} \\
  chr1:100-200 \\
  control-statistics.zip \\
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
"""Compute summary statistics for control gene from BAM
files."""
    )
    parser.add_argument(
        'gene',
        help=
"""Control gene (recommended choices: 'EGFR', 'RYR1',
'VDR'). Alternatively, you can provide a custom region
(format: chrom:start-end)."""
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
to indicate probed regions."""
    )

def main(args):
    result = utils.compute_control_statistics(
        args.gene, args.bams, assembly=args.assembly, bed=args.bed
    )
    result.to_file(args.control_statistics)
