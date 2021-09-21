import sys

from ..api import pipeline

import fuc

description = f"""
#########################################
# Run NGS pipeline for the target gene. #
#########################################

Usage examples:
  $ pypgx {fuc.api.common._script_name()} CYP2D6 CYP2D6-pipeline --vcf input.vcf --panel ref.vcf --tsv input.tsv --control-statistcs control-statistics-VDR.zip
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
        '--vcf',
        metavar='PATH',
        help='VCF file.'
    )
    parser.add_argument(
        '--panel',
        metavar='PATH',
        help='Reference haplotype panel.'
    )
    parser.add_argument(
        '--tsv',
        metavar='PATH',
        help='TSV file containing read depth (zipped or unzipped).'
    )
    parser.add_argument(
        '--control-statistics',
        metavar='PATH',
        help='Archive file with the semandtic type SampleTable[Statistcs].'
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
        args.gene, args.output, vcf=args.vcf, panel=args.panel, tsv=args.tsv,
        control_statistics=args.control_statistics, force=args.force,
        do_not_plot_copy_number=args.do_not_plot_copy_number,
        do_not_plot_allele_fraction=args.do_not_plot_allele_fraction
    )
