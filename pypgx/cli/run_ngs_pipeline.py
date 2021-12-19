import sys

from ..api import pipeline

import fuc

description = """
Run PyPGx's genotyping pipeline for NGS data.

During copy number analysis, if the input data is targeted sequencing, the
command will apply inter-sample normalization using summary statistics across
all samples. For best results, it is recommended to specify known samples
without SV using --samples-without-sv.
"""

epilog = f"""
[Example] To genotype the CYP3A5 gene, which does not have SV, from WGS data:
  $ pypgx {fuc.api.common._script_name()} \\
  CYP3A5 \\
  CYP3A5-pipeline \\
  --variants variants.vcf.gz

[Example] To genotype the CYP2D6 gene, which does have SV, from WGS data:
  $ pypgx {fuc.api.common._script_name()} \\
  CYP2D6 \\
  CYP2D6-pipeline \\
  --variants variants.vcf.gz \\
  --depth-of-coverage depth-of-coverage.tsv \\
  --control-statistcs control-statistics-VDR.zip

[Example] To genotype the CYP2D6 gene from targeted sequencing data:
  $ pypgx {fuc.api.common._script_name()} \\
  CYP2D6 \\
  CYP2D6-pipeline \\
  --variants variants.vcf.gz \\
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
        help='Input VCF file must be already BGZF compressed (.gz) \n'
             'and indexed (.tbi) to allow random access. Statistical \n'
             'haplotype phasing will be skipped if input VCF is \n'
             'already fully phased.'
    )
    parser.add_argument(
        '--depth-of-coverage',
        metavar='PATH',
        help='Depth of coverage file (compressed or uncompressed).'
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
        help="Reference genome assembly (default: 'GRCh37') \n"
             "(choices: 'GRCh37', 'GRCh38')."
    )
    parser.add_argument(
        '--panel',
        metavar='PATH',
        help='VCF file corresponding to a reference haplotype panel \n'
             '(compressed or uncompressed). By default, the 1KGP \n'
             'panel is used.'
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
        help='Specify which samples should be included for analysis \n'
             'by providing a text file (.txt, .tsv, .csv, or .list) \n'
             'containing one sample per line. Alternatively, you can \n'
             'provide a list of samples.'
    )
    parser.add_argument(
        '--exclude',
        action='store_true',
        help='Exclude specified samples.'
    )
    parser.add_argument(
        '--samples-without-sv',
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
        exclude=args.exclude, samples_without_sv=args.samples_without_sv,
        do_not_plot_copy_number=args.do_not_plot_copy_number,
        do_not_plot_allele_fraction=args.do_not_plot_allele_fraction,
        platform=args.platform
    )
