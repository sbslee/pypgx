import argparse
import pypgx
from .cli import commands
from .version import __version__

def _func(l):
    """Convert a list to a pretty text."""
    return '{' + ", ".join([f"'{x}'" for x in l]) + '}'

def _get_parser():
    parser = argparse.ArgumentParser(add_help=False)

    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="Show the version and exit."
    )

    parser.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit."
    )

    subparsers = parser.add_subparsers(
        dest="command",
        metavar="COMMAND",
        help="Name of the command."
    )

    subparsers.required = True

    compare_stargazer_calls_parser = subparsers.add_parser(
        "compare-stargazer-calls",
        add_help=False,
        help=("Compute the concordance between two 'genotype-calls.tsv' "
              "files created by Stargazer."),
        description=("Compute the concordance between two "
                     "'genotype-calls.tsv' files created by Stargazer.")
    )

    compare_stargazer_calls_parser._optionals.title = "Arguments"

    compare_stargazer_calls_parser.add_argument(
        "-r",
        "--ref-file",
        metavar="PATH",
        required=True,
        help=("Path to the reference or truth 'genotype-calls.tsv' file "
              "created by Stargazer. [required]")
    )

    compare_stargazer_calls_parser.add_argument(
        "-t",
        "--test-file",
        metavar="PATH",
        required=True,
        help=("Path to the test 'genotype-calls.tsv' file "
              "created by Stargazer. [required]")
    )

    compare_stargazer_calls_parser.add_argument(
        "-o",
        "--output-file",
        metavar="PATH",
        required=True,
        help=("Path to the output file. [required]")
    )

    compare_stargazer_calls_parser.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit."
    )

    calculate_read_depth_parser = subparsers.add_parser(
        "calculate-read-depth",
        add_help=False,
        help=("Create a GDF (GATK DepthOfCoverage Format) file for "
              "Stargazer from BAM files by computing read depth."),
        description=("Create a GDF (GATK DepthOfCoverage Format) file for "
                     "Stargazer from BAM files by computing read depth.")
    )

    calculate_read_depth_parser._optionals.title = "Arguments"

    calculate_read_depth_parser.add_argument(
        "-t",
        "--target-gene",
        metavar="TEXT",
        required=True,
        help=("Name of the target gene. Choices: "
              f"{_func(pypgx.target_genes)}. [required]")
    )

    calculate_read_depth_parser.add_argument(
        "-c",
        "--control-gene",
        metavar="TEXT",
        required=True,
        help=("Name of a preselected control gene. Used for intrasample "
              "normalization during copy number analysis by Stargazer. "
              f"Choices: {_func(pypgx.control_genes)}. Alternatively, "
              "you can provide a custom genomic region with the "
              "'chr:start-end' format (e.g. chr12:48232319-48301814). "
              "[required]")
    )

    calculate_read_depth_parser.add_argument(
        "-i",
        "--bam-path",
        metavar="PATH",
        required=True,
        help=("Read BAM files from PATH, one file path per line. "
              "Also accepts single BAM file. [required]")
    )

    calculate_read_depth_parser.add_argument(
        "-o",
        "--output-file",
        metavar="PATH",
        required=True,
        help=("Path to the output file. [required]")
    )

    calculate_read_depth_parser.add_argument(
        "-a",
        "--genome-build",
        metavar="TEXT",
        default="hg19",
        help=("Build of the reference genome assembly. "
              "Choices: {'hg19', 'hg38'}. [default: 'hg19']")
    )

    calculate_read_depth_parser.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit."
    )

    call_variants_gatk_sge_parser = subparsers.add_parser(
        "call-variants-gatk-sge",
        add_help=False,
        help=("Create a VCF (Variant Call Format) file for Stargazer "
              "from BAM files by calling SNVs and indels."),
        description=("Create a VCF (Variant Call Format) file for Stargazer "
                     "from BAM files by calling SNVs and indels.")
    )

    call_variants_gatk_sge_parser._optionals.title = "Arguments"

    call_variants_gatk_sge_parser.add_argument(
        "-t",
        "--target-gene",
        metavar="TEXT",
        required=True,
        help=("Name of the target gene. Choices: "
              f"{_func(pypgx.target_genes)}. [required]")
    )

    call_variants_gatk_sge_parser.add_argument(
        "-i",
        "--bam-path",
        metavar="PATH",
        required=True,
        help=("Read BAM files from PATH, one file path per line. [required]")
    )

    call_variants_gatk_sge_parser.add_argument(
        "-f",
        "--fasta-file",
        metavar="PATH",
        required=True,
        help=("Path to a reference FASTA file. [required]")
    )

    call_variants_gatk_sge_parser.add_argument(
        "-o",
        "--output-dir",
        metavar="PATH",
        required=True,
        help=("Path to the output directory. [required]")
    )

    call_variants_gatk_sge_parser.add_argument(
        "-a",
        "--genome-build",
        metavar="TEXT",
        default="hg19",
        help=("Build of the reference genome assembly. "
              "Choices: {'hg19', 'hg38'}. [default: 'hg19']")
    )

    call_variants_gatk_sge_parser.add_argument(
        "-d",
        "--dbsnp-file",
        metavar="PATH",
        help=("Path to a dbSNP file (.vcf or .vcf.gz). Used to assign "
              "rs ID to observed variants.")
    )

    call_variants_gatk_sge_parser.add_argument(
        "-j",
        "--java-options",
        metavar="TEXT",
        help=("Options passed to Java to run GATK. Must be a quoted string "
              'proceeded by an equal sign (e.g. -j="-Xmx4G").')
    )

    call_variants_gatk_sge_parser.add_argument(
        "-q",
        "--qsub-options",
        metavar="TEXT",
        help=("Options passed to SGE. Must be a quoted string proceeded "
              'by an equal sign (e.g. -q="-l mem_requested=4G").')
    )

    call_variants_gatk_sge_parser.add_argument(
        "-c",
        "--conda-env",
        metavar="TEXT",
        help=("Name of the conda environment to be activated when the jobs "
              "are submitted to SGE.")
    )

    call_variants_gatk_sge_parser.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit."
    )

    return parser

def main():
    parser = _get_parser()
    args = parser.parse_args()
    command = args.command
    delattr(args, "command")
    commands[command](**vars(args))

if __name__ == "__main__":
    main()
