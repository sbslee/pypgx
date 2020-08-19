import sys
import argparse
import timeit
import datetime

from .version import __version__
from .common import (
    get_logger,
    get_file_list,
    read_file_list,
)

from .bam2gt import bam2gt
from .bam2gt2 import bam2gt2
from .gt2pt import gt2pt
from .bam2vcf import bam2vcf
from .bam2vcf2 import bam2vcf2
from .bam2gdf import bam2gdf
from .gt2html import gt2html
from .bam2html import bam2html
from .fq2bam import fq2bam
from .bam2bam import bam2bam
from .bam2sdf import bam2sdf
from .sdf2gdf import sdf2gdf
from .pgkb import pgkb
from .minivcf import minivcf
from .mergevcf import mergevcf
from .summary import summary
from .meta import meta
from .compare import compare
from .check import check
from .liftover import liftover
from .peek import peek
from .viewsnp import viewsnp
from .compgt import compgt
from .compvcf import compvcf
from .unicov import unicov

PYPGX_TOOLS = {
    "bam2gt": bam2gt,
    "bam2gt2": bam2gt2,
    "gt2pt": gt2pt,
    "bam2vcf": bam2vcf,
    "bam2vcf2": bam2vcf2,
    "bam2gdf": bam2gdf,
    "gt2html": gt2html,
    "bam2html": bam2html,
    "fq2bam": fq2bam,
    "bam2bam": bam2bam,
    "bam2sdf": bam2sdf,
    "sdf2gdf": sdf2gdf,
    "pgkb": pgkb,
    "minivcf": minivcf,
    "mergevcf": mergevcf,
    "summary": summary,
    "meta": meta,
    "compare": compare,
    "check": check,
    "liftover": liftover,
    "peek": peek,
    "viewsnp": viewsnp,
    "compgt": compgt,
    "compvcf": compvcf,
    "unicov": unicov,
}

def get_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"PyPGx v{__version__}",
        help="print the PyPGx version number and exit"
    )

    output_parser = argparse.ArgumentParser(add_help=False)
    output_parser.add_argument(
        "-o",
        "--output",
        metavar="FILE",
        help="write output to FILE [stdout]"
    )

    bam_getter_parser = argparse.ArgumentParser(add_help=False)
    bam_getter_parser.add_argument(
        "--bam_dir",
        metavar="DIR",
        help="treat any BAM files in DIR as input"
    )
    bam_getter_parser.add_argument(
        "--bam_list",
        metavar="FILE",
        help="read BAM files from FILE, one file path per line"
    )

    subparsers = parser.add_subparsers(
        dest="tool",
        metavar="tool",
        help="name of the tool",
    )

    bam2gt_parser = subparsers.add_parser(
        "bam2gt",
        help="convert BAM files to a genotype file",
        parents=[bam_getter_parser]
    )
    bam2gt_parser.add_argument(
        "snp_caller",
        help="SNP caller ('gatk' or 'bcftools')"
    )
    bam2gt_parser.add_argument(
        "fasta_file",
        help="reference FASTA file"
    )
    bam2gt_parser.add_argument(
        "target_gene",
        help="name of target gene (e.g. 'cyp2d6')"
    )
    bam2gt_parser.add_argument(
        "genome_build",
        help="genome build ('hg19' or 'hg38')"
    )
    bam2gt_parser.add_argument(
        "data_type",
        help="type of sequencing data ('wgs' or 'ts')",
    )
    bam2gt_parser.add_argument(
        "proj_dir",
        help="output files will be written to this directory",
    )
    bam2gt_parser.add_argument(
        "bam_file",
        nargs="*",
        help="input BAM files"
    )
    bam2gt_parser.add_argument(
        "--control_gene",
        metavar="STR",
        help="name or region of control gene (e.g. ‘vdr’, "
            + "‘chr12:48232319-48301814’)"
    )
    bam2gt_parser.add_argument(
        "--dbsnp_file",
        metavar="FILE",
        help="dbSNP VCF file, used by GATK to add rs numbers"
    )
    bam2gt_parser.add_argument(
        "--temp_dir",
        metavar="DIR",
        help="temporary files will be written to this directory"
    )
    bam2gt_parser.add_argument(
        "--plot",
        action="store_true",
        help="output copy number plots",
    )

    bam2gt2_parser = subparsers.add_parser(
        "bam2gt2",
        help="convert BAM files to genotype files [SGE]",
        parents=[output_parser]
    )
    bam2gt2_parser.add_argument(
        "conf_file",
        help="configuration file",
    )

    gt2pt_parser = subparsers.add_parser(
        "gt2pt",
        help="convert a genotype file to phenotypes",
        parents=[output_parser]
    )
    gt2pt_parser.add_argument(
        "gt_file",
        help="genotype file from Stargazer ('genotype.txt')",
    )

    bam2vcf_parser = subparsers.add_parser(
        "bam2vcf",
        help="convert BAM files to a VCF file",
        parents=[output_parser]
    )
    bam2vcf_parser.add_argument(
        "snp_caller",
        help="SNP caller ('gatk' or 'bcftools')"
    )
    bam2vcf_parser.add_argument(
        "fasta_file",
        help="reference FASTA file"
    )
    bam2vcf_parser.add_argument(
        "target_gene",
        help="name or region of target gene (e.g. ‘cyp2d6’, "
            + "‘chr22:42512500-42551883’)"
    )
    bam2vcf_parser.add_argument(
        "output_file",
        help="write output to this file"
    )
    bam2vcf_parser.add_argument(
        "genome_build",
        help="genome build ('hg19' or 'hg38')"
    )
    bam2vcf_parser.add_argument(
        "bam_file",
        nargs="*",
        help="input BAM files"
    )
    bam2vcf_parser.add_argument(
        "--bam_dir",
        metavar="DIR",
        help="use all BAM files in this directory as input"
    )
    bam2vcf_parser.add_argument(
        "--bam_list",
        metavar="FILE",
        help="list of input BAM files, one file per line"
    )
    bam2vcf_parser.add_argument(
        "--dbsnp_file",
        metavar="FILE",
        help="dbSNP VCF file, used by GATK to add rs numbers"
    )
    bam2vcf_parser.add_argument(
        "--java_options",
        metavar="STR",
        help="Java-specific arguments for GATK (e.g. '-Xmx4G')"
    )
    bam2vcf_parser.add_argument(
        "--temp_dir",
        metavar="DIR",
        help="temporary files will be written to this directory"
    )

    bam2vcf2_parser = subparsers.add_parser(
        "bam2vcf2",
        help="convert BAM files to a VCF file [SGE]",
        parents=[output_parser]
    )
    bam2vcf2_parser.add_argument(
        "conf_file",
        help="configuration file",
    )

    bam2gdf_parser = subparsers.add_parser(
        "bam2gdf",
        help="convert BAM files to a GDF file",
        parents=[output_parser]
    )
    bam2gdf_parser.add_argument(
        "genome_build",
        help="genome build ('hg19' or 'hg38')",
    )
    bam2gdf_parser.add_argument(
        "target_gene",
        help="name of target gene (e.g. 'cyp2d6')",
    )
    bam2gdf_parser.add_argument(
        "control_gene",
        help="name or region of control gene (e.g. ‘vdr’, "
            + "‘chr12:48232319-48301814’)"
    )
    bam2gdf_parser.add_argument(
        "output_file",
        help="write output to this file"
    )
    bam2gdf_parser.add_argument(
        "bam_file",
        nargs="*",
        help="input BAM files"
    )
    bam2gdf_parser.add_argument(
        "--bam_dir",
        metavar="DIR",
        help="use all BAM files in this directory as input"
    )
    bam2gdf_parser.add_argument(
        "--bam_list",
        metavar="FILE",
        help="list of input BAM files, one file per line"
    )

    gt2html_parser = subparsers.add_parser(
        "gt2html",
        help="convert a genotype file to an HTML report",
        parents=[output_parser]
    )
    gt2html_parser.add_argument(
        "gt_file",
        help="genotype file from Stargazer ('genotype.txt')",
    )

    bam2html_parser = subparsers.add_parser(
        "bam2html",
        help="convert a BAM file to an HTML report [SGE]",
        parents=[output_parser]
    )
    bam2html_parser.add_argument(
        "conf_file",
        help="configuration file",
    )

    fq2bam_parser = subparsers.add_parser(
        "fq2bam",
        help="convert FASTQ files to BAM files [SGE]",
        parents=[output_parser]
    )
    fq2bam_parser.add_argument(
        "conf_file",
        help="configuration file",
    )

    bam2bam_parser = subparsers.add_parser(
        "bam2bam",
        help="realign BAM files to another reference genome [SGE]",
        parents=[output_parser]
    )
    bam2bam_parser.add_argument(
        "conf_file",
        help="configuration file",
    )

    bam2sdf_parser = subparsers.add_parser(
        "bam2sdf",
        help="convert BAM files to a SDF file",
        parents=[output_parser]
    )
    bam2sdf_parser.add_argument(
        "gb",
        help="genome build",
    )
    bam2sdf_parser.add_argument(
        "tg",
        help="target gene",
    )
    bam2sdf_parser.add_argument(
        "cg",
        help="control gene or region",
    )
    bam2sdf_parser.add_argument(
        "bam",
        nargs="+",
        help="BAM file",
    )

    sdf2gdf_parser = subparsers.add_parser(
        "sdf2gdf",
        help="convert a SDF file to a GDF file",
        parents=[output_parser]
    )
    sdf2gdf_parser.add_argument(
        "sdf",
        help="SDF file",
    )
    sdf2gdf_parser.add_argument(
        "id",
        nargs="+",
        help="sample ID",
    )

    pgkb_parser = subparsers.add_parser(
        "pgkb",
        help="extract CPIC guidelines using PharmGKB API",
        parents=[output_parser]
    )
    pgkb_parser.add_argument(
        "--test_mode",
        action="store_true",
        help="only extract first three guidelines for testing",
    )

    minivcf_parser = subparsers.add_parser(
        "minivcf",
        help="slice VCF file",
        parents=[output_parser]
    )
    minivcf_parser.add_argument(
        "vcf_file",
        help="VCF file",
    )
    minivcf_parser.add_argument(
        "region",
        help="target region",
    )

    mergevcf_parser = subparsers.add_parser(
        "mergevcf",
        help="merge VCF files",
        parents=[output_parser]
    )
    mergevcf_parser.add_argument(
        "vcf_file",
        nargs="+",
        help="VCF files to be merged",
    )
    mergevcf_parser.add_argument(
        "--region",
        metavar="STR",
        help="target region",
    )

    summary_parser = subparsers.add_parser(
        "summary",
        help="create summary file using Stargazer data",
        parents=[output_parser]
    )
    summary_parser.add_argument(
        "gt",
        help="genotype file",
    )

    meta_parser = subparsers.add_parser(
        "meta",
        help="create meta file from summary files",
        parents=[output_parser]
    )
    meta_parser.add_argument(
        "sf",
        nargs="+",
        help="summary file",
    )

    compare_parser = subparsers.add_parser(
        "compare",
        help="compare genotype files",
        parents=[output_parser]
    )
    compare_parser.add_argument(
        "gt",
        nargs="+",
        help="genotype file",
    )

    check_parser = subparsers.add_parser(
        "check",
        help="check table files for Stargazer",
        parents=[output_parser]
    )
    check_parser.add_argument(
        "star_table",
        help="star allele table file",
    )
    check_parser.add_argument(
        "snp_table",
        help="SNP table file",
    )

    liftover_parser = subparsers.add_parser(
        "liftover",
        help="convert variants in SNP table from hg19 to hg38",
        parents=[output_parser]
    )
    liftover_parser.add_argument(
        "star_table",
        help="star allele table file",
    )
    liftover_parser.add_argument(
        "snp_table",
        help="SNP table file",
    )
    liftover_parser.add_argument(
        "target_gene",
        help="target gene",
    )

    peek_parser = subparsers.add_parser(
        "peek",
        help="find all possible star alleles from VCF file",
        parents=[output_parser]
    )
    peek_parser.add_argument(
        "vcf_file",
        help="Stargazer VCF file (finalized.vcf)",
    )

    viewsnp_parser = subparsers.add_parser(
        "viewsnp",
        help="view SNP data for pairs of sample/star allele",
        parents=[output_parser]
    )
    viewsnp_parser.add_argument(
        "vcf_file",
        help="Stargazer VCF file ('finalized.vcf')",
    )
    viewsnp_parser.add_argument(
        "query",
        nargs="+",
        help="pair of sample and star allele separated by '/' "
            + "(e.g. 'SAMPLE1/*4')",
    )

    compgt_parser = subparsers.add_parser(
        "compgt",
        help="compute the concordance between two genotype files",
        parents=[output_parser]
    )
    compgt_parser.add_argument(
        "truth_file",
        help="truth genotype file",
    )
    compgt_parser.add_argument(
        "test_file",
        help="test genotype file",
    )
    compgt_parser.add_argument(
        "sample_map",
        help="tab-delimited text file with two columns representing "
            + "the truth and test sample names",
    )

    compvcf_parser = subparsers.add_parser(
        "compvcf",
        help="compute the concordance between two VCF files",
        parents=[output_parser]
    )
    compvcf_parser.add_argument(
        "truth_file",
        help="truth VCF file",
    )
    compvcf_parser.add_argument(
        "test_file",
        help="test VCF file",
    )
    compvcf_parser.add_argument(
        "sample_map",
        help="tab-delimited text file with two columns representing "
            + "the truth and test sample names",
    )

    unicov_parser = subparsers.add_parser(
        "unicov",
        help="compute the uniformity of sequencing coverage",
        parents=[output_parser]
    )
    unicov_parser.add_argument(
        "bed_file",
        help="BED file",
    )
    unicov_parser.add_argument(
        "bam_file",
        nargs="*",
        help="input BAM files"
    )
    unicov_parser.add_argument(
        "--bam_dir",
        metavar="DIR",
        help="use all BAM files in this directory as input"
    )
    unicov_parser.add_argument(
        "--bam_list",
        metavar="FILE",
        help="list of input BAM files, one file per line"
    )

    return parser

def main():
    start_time = timeit.default_timer()
    parser = get_parser()
    args = parser.parse_args()

    logger = get_logger()
    logger.info(f"PyPGx v{__version__}")
    logger.info("Command:")
    logger.info("    {}".format(" ".join(sys.argv)))

    result = PYPGX_TOOLS[args.tool](**vars(args))

    stop_time = timeit.default_timer()
    elapsed_time = str(
        datetime.timedelta(seconds=(stop_time - start_time))
    ).split(".")[0]

    logger.info(f"Elapsed time: {elapsed_time}")
    logger.info("PyPGx finished")

    if result:
        if args.output:
            with open(args.output, "w") as f:
                f.write(result)
        else:
            sys.stdout.write(result)

if __name__ == "__main__":
    main()
