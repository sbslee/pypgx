import sys

from ..api import pipeline

import fuc

description = """
Run PyPGx's genotyping pipeline for NGS data.

During copy number analysis, if the input data is targeted sequencing, the
command will apply inter-sample normalization using summary statistics across
all samples. For best results, it is recommended to specify known samples
without SV using --samples.
"""

epilog = f"""
[Example] To genotype the CYP3A5 gene, which does not have SV, from WGS data:
  $ pypgx {fuc.api.common._script_name()} \\
  CYP3A5 \\
  CYP3A5-pipeline \\
  --variants variants.vcf

[Example] To genotype the CYP2D6 gene, which does have SV, from WGS data:
  $ pypgx {fuc.api.common._script_name()} \\
  CYP2D6 \\
  CYP2D6-pipeline \\
  --variants variants.vcf \\
  --depth-of-coverage depth-of-coverage.tsv \\
  --control-statistcs control-statistics-VDR.zip

[Example] To genotype the CYP2D6 gene from targeted sequencing data:
  $ pypgx {fuc.api.common._script_name()} \\
  CYP2D6 \\
  CYP2D6-pipeline \\
  --variants variants.vcf \\
  --depth-of-coverage depth-of-coverage.tsv \\
  --control-statistcs control-statistics-VDR.zip \\
  --platform Targeted
"""

def create_parser(subparsers):
    parser = fuc.api.common._add_parser(
        subparsers,
        fuc.api.common._script_name(),
        description=description,
        epilog=epilog,
        help="Run PyPGx's genotyping pipeline for NGS data.",
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
        help='Archive file with the semandtic type \n'
             'SampleTable[Statistcs].'
    )
    parser.add_argument(
        '--platform',
        metavar='TEXT',
        default='WGS',
        choices=['WGS', 'Targeted'],
        help="Genotyping platform (default: 'WGS') (choices: 'WGS', \n"
             "'Targeted')"
    )
    parser.add_argument(
        '--assembly',
        metavar='TEXT',
        default='GRCh37',
        help="Reference genome assembly (default: 'GRCh37') (choices: \n"
             "'GRCh37', 'GRCh38')."
    )
    parser.add_argument(
        '--panel',
        metavar='PATH',
        help='VCF file corresponding to a reference haplotype panel \n'
             '(zipped or unzipped). By default, the 1KGP panel is \n'
             'used.'
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help='Overwrite output directory if it already exists.'
    )
    parser.add_argument(
        '--samples',
        metavar='TEXT',
        nargs='+',
        help="List of known samples without SV."
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
        control_statistics=args.control_statistics, assembly=args.assembly, 
        panel=args.panel, force=args.force, samples=args.samples,
        do_not_plot_copy_number=args.do_not_plot_copy_number,
        do_not_plot_allele_fraction=args.do_not_plot_allele_fraction,
        platform=args.platform
    )
