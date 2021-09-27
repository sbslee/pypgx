import sys

from ..api import pipeline

import fuc

description = f"""
#########################################
# Run NGS pipeline for the target gene. #
#########################################

Usage examples:
  $ pypgx {fuc.api.common._script_name()} CYP2D6 CYP2D6-pipeline --variants variants.vcf --depth-of-coverage depth-of-coverage.tsv --control-statistcs control-statistics-VDR.zip
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        help='Run NGS pipeline for the target gene.',
        description=description,
    )
    parser.add_argument(
        'gene',
        help='Target gene.'
    )
    parser.add_argument(
        'output',
        help='Output directory.'
    )
    parser.add_argument(
        '--variants',
        metavar='PATH',
        help='VCF file (zipped or unzipped).'
    )
    parser.add_argument(
        '--depth-of-coverage',
        metavar='PATH',
        help='Depth of coverage file (zipped or unzipped).'
    )
    parser.add_argument(
        '--control-statistics',
        metavar='PATH',
        help='Archive file with the semandtic type SampleTable[Statistcs].'
    )
    parser.add_argument(
        '--panel',
        metavar='PATH',
        help='Reference haplotype panel. By default, the 1KGP panel is used.'
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help='Overwrite output directory if it already exists.'
    )
    parser.add_argument(
        '--do-not-plot-copy-number',
        action='store_true',
        help='Do not plot copy number profile.'
    )
    parser.add_argument(
        '--do-not-plot-allele-fraction',
        action='store_true',
        help='Do not plot allele fraction profile.'
    )

def main(args):
    pipeline.run_ngs_pipeline(
        args.gene, args.output, variants=args.variants,
        depth_of_coverage=args.depth_of_coverage,
        control_statistics=args.control_statistics,
        panel=args.panel, force=args.force,
        do_not_plot_copy_number=args.do_not_plot_copy_number,
        do_not_plot_allele_fraction=args.do_not_plot_allele_fraction
    )
