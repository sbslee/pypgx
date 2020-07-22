import sys
import argparse
import timeit
import datetime

from .version import __version__
from .common import (
    logging,
    get_file_list,
    read_file_list,
)

from .pgkb import pgkb
from .report import report
from .sdf2gdf import sdf2gdf
from .bam2sdf import bam2sdf
from .bam2gdf import bam2gdf
from .minivcf import minivcf
from .merge import merge
from .summary import summary
from .meta import meta
from .compare import compare
from .remap import remap
from .fq2bam import fq2bam
from .sges import sges
from .sgep import sgep
from .sgea import sgea
from .cpa import cpa
from .plotcov import plotcov
from .check import check
from .liftover import liftover
from .peek import peek
from .snp import snp
from .bam2vcf import bam2vcf
from .genotype import genotype
from .gt2pt import gt2pt

def _get_bam_list(bc, bd, bl):
    """Get the list of BAM files.

    Returns:
        list[str]: List of BAM files.

    Args:
        bc (list[str]): List of BAM files from the console.
        bd (str): Path to directory containing BAM files.
        bl (str): Path to file containing BAM filenames.
    """
    if bc:
        result = bc
    elif bd:
        result = get_file_list(bd, ".bam")
    elif bl:
        result = read_file_list(bl)
    else:
        raise ValueError("BAM files not provided")
    return result

def get_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--version",
        action="version",
        version=f"PyPGx v{__version__}",
        help="print the PyPGx version number and exit"
    )

    subparsers = parser.add_subparsers(
        dest="tool",
        metavar="tool",
        help="name of the tool",
    )

    pgkb_parser = subparsers.add_parser(
        "pgkb",
        help="extract CPIC guidelines using PharmGKB API",
    )
    pgkb_parser.add_argument(
        "-o",
        metavar="FILE",
        help="output to FILE [stdout]",
    )
    pgkb_parser.add_argument(
        "-t",
        action="store_true",
        help="extract first three guidelines for testing",
    )

    report_parser = subparsers.add_parser(
        "report",
        help="create HTML report using Stargazer data",
    )
    report_parser.add_argument(
        "gt",
        help="genotype file",
    )
    report_parser.add_argument(
        "-o",
        metavar="FILE",
        help="output to FILE [stdout]",
    )

    sdf2gdf_parser = subparsers.add_parser(
        "sdf2gdf",
        help="create GDF file from SDF file",
    )
    sdf2gdf_parser.add_argument(
        "sdf",
        help="SDF file",
    )
    sdf2gdf_parser.add_argument(
        "-o",
        metavar="FILE",
        help="output to FILE [stdout]",
    )
    sdf2gdf_parser.add_argument(
        "id",
        nargs="+",
        help="sample ID",
    )

    bam2sdf_parser = subparsers.add_parser(
        "bam2sdf",
        help="create SDF file from BAM file(s)",
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
    bam2sdf_parser.add_argument(
        "-o",
        metavar="FILE",
        help="output to FILE [stdout]",
    )

    bam2gdf_parser = subparsers.add_parser(
        "bam2gdf",
        help="create GDF file from BAM file(s)",
    )
    bam2gdf_parser.add_argument(
        "gb",
        help="genome build",
    )
    bam2gdf_parser.add_argument(
        "tg",
        help="target gene",
    )
    bam2gdf_parser.add_argument(
        "cg",
        help="control gene or region",
    )
    bam2gdf_parser.add_argument(
        "bam",
        nargs="*",
        help="BAM file",
    )
    bam2gdf_parser.add_argument(
        "--bd",
        metavar="DIR",
        help="directory containing BAM files",
    )
    bam2gdf_parser.add_argument(
        "--bl",
        metavar="FILE",
        help="list of BAM files, one file per line"
    )
    bam2gdf_parser.add_argument(
        "-o",
        metavar="FILE",
        help="output to FILE [stdout]",
    )

    minivcf_parser = subparsers.add_parser(
        "minivcf",
        help="slice VCF file",
    )
    minivcf_parser.add_argument(
        "vcf",
        help="VCF file",
    )
    minivcf_parser.add_argument(
        "region",
        help="target region",
    )
    minivcf_parser.add_argument(
        "-o",
        metavar="FILE",
        help="output to FILE [stdout]",
    )

    merge_parser = subparsers.add_parser(
        "merge",
        help="merge VCF files",
    )
    merge_parser.add_argument(
        "vcf",
        nargs="+",
        help="VCF file",
    )
    merge_parser.add_argument(
        "-r",
        metavar="STR",
        help="target region",
    )
    merge_parser.add_argument(
        "-o",
        metavar="FILE",
        help="output to FILE [stdout]",
    )

    summary_parser = subparsers.add_parser(
        "summary",
        help="create summary file using Stargazer data",
    )
    summary_parser.add_argument(
        "gt",
        help="genotype file",
    )
    summary_parser.add_argument(
        "-o",
        metavar="FILE",
        help="output to FILE [stdout]",
    )

    meta_parser = subparsers.add_parser(
        "meta",
        help="create meta file from summary files",
    )
    meta_parser.add_argument(
        "sf",
        nargs="+",
        help="summary file",
    )
    meta_parser.add_argument(
        "-o",
        metavar="FILE",
        help="output to FILE [stdout]",
    )

    compare_parser = subparsers.add_parser(
        "compare",
        help="compare genotype files",
    )
    compare_parser.add_argument(
        "gt",
        nargs="+",
        help="genotype file",
    )
    compare_parser.add_argument(
        "-o",
        metavar="FILE",
        help="output to FILE [stdout]",
    )

    remap_parser = subparsers.add_parser(
        "remap",
        help="remap BAM file(s) to different reference",
    )
    remap_parser.add_argument(
        "conf",
        help="configuration file",
    )

    fq2bam_parser = subparsers.add_parser(
        "fq2bam",
        help="create BAM file(s) from FASTQ file(s)",
    )
    fq2bam_parser.add_argument(
        "conf",
        help="configuration file",
    )

    sges_parser = subparsers.add_parser(
        "sges",
        help="run per-sample genotyping with Stargazer",
    )
    sges_parser.add_argument(
        "conf",
        help="configuration file",
    )

    sgep_parser = subparsers.add_parser(
        "sgep",
        help="run per-project genotyping with Stargazer (1)",
    )
    sgep_parser.add_argument(
        "conf",
        help="configuration file",
    )

    sgea_parser = subparsers.add_parser(
        "sgea",
        help="run per-project genotyping with Stargazer (2)",
    )
    sgea_parser.add_argument(
        "conf",
        help="configuration file",
    )

    cpa_parser = subparsers.add_parser(
        "cpa",
        help="run change point analysis for copy number",
    )
    cpa_parser.add_argument(
        "rdata",
        help="RData file",
    )
    cpa_parser.add_argument(
        "-o",
        metavar="FILE",
        help="output to FILE [stdout]",
    )

    plotcov_parser = subparsers.add_parser(
        "plotcov",
        help="plot coverage data to PDF file",
    )
    plotcov_parser.add_argument(
        "sdf",
        help="SDF file",
    )
    plotcov_parser.add_argument(
        "out",
        help="PDF file",
    )

    check_parser = subparsers.add_parser(
        "check",
        help="check table files for Stargazer",
    )
    check_parser.add_argument(
        "star",
        help="star allele table file",
    )
    check_parser.add_argument(
        "snp",
        help="SNP table file",
    )

    liftover_parser = subparsers.add_parser(
        "liftover",
        help="convert variants in SNP table from hg19 to hg38",
    )
    liftover_parser.add_argument(
        "star",
        help="star allele table file",
    )
    liftover_parser.add_argument(
        "snp",
        help="SNP table file",
    )
    liftover_parser.add_argument(
        "tg",
        help="target gene",
    )
    liftover_parser.add_argument(
        "-o",
        metavar="FILE",
        help="output to FILE [stdout]",
    )

    peek_parser = subparsers.add_parser(
        "peek",
        help="find all possible star alleles from VCF file",
    )
    peek_parser.add_argument(
        "vcf",
        help="VCF file",
    )
    peek_parser.add_argument(
        "-o",
        metavar="FILE",
        help="output to FILE [stdout]",
    )

    snp_parser = subparsers.add_parser(
        "snp",
        help="view variant data for sample/star allele pairs",
    )
    snp_parser.add_argument(
        "vcf",
        help="VCF file",
    )
    snp_parser.add_argument(
        "pair",
        nargs="+",
        help="sample/star allele pair",
    )
    snp_parser.add_argument(
        "-o",
        metavar="FILE",
        help="output to FILE [stdout]",
    )

    bam2vcf_parser = subparsers.add_parser(
        "bam2vcf",
        help="create VCF file from BAM file(s)",
    )
    bam2vcf_parser.add_argument(
        "gb",
        help="genome build",
    )
    bam2vcf_parser.add_argument(
        "tg",
        help="target gene or region"
    )
    bam2vcf_parser.add_argument(
        "fa",
        help="FASTA file"
    )
    bam2vcf_parser.add_argument(
        "bam",
        nargs="*",
        help="BAM file",
    )
    bam2vcf_parser.add_argument(
        "--bd",
        metavar="DIR",
        help="directory containing BAM files",
    )
    bam2vcf_parser.add_argument(
        "--bl",
        metavar="FILE",
        help="list of BAM files, one file per line"
    )
    bam2vcf_parser.add_argument(
        "-o",
        metavar="FILE",
        help="output to FILE [stdout]",
    )

    genotype_parser = subparsers.add_parser(
        "genotype",
        help="call star alleles from BAM file(s)",
    )
    genotype_parser.add_argument(
        "fa",
        help="FASTA file"
    )
    genotype_parser.add_argument(
        "dt",
        help="input data type",
    )
    genotype_parser.add_argument(
        "gb",
        help="genome build",
    )
    genotype_parser.add_argument(
        "tg",
        help="target gene",
    )
    genotype_parser.add_argument(
        "out",
        help="output project directory",
    )
    genotype_parser.add_argument(
        "bam",
        nargs="*",
        help="BAM file",
    )
    genotype_parser.add_argument(
        "--bd",
        metavar="DIR",
        help="directory containing BAM files",
    )
    genotype_parser.add_argument(
        "--bl",
        metavar="FILE",
        help="list of BAM files, one file per line"
    )
    genotype_parser.add_argument(
        "--cg",
        metavar="STR",
        help="control gene",
    )

    gt2pt_parser = subparsers.add_parser(
        "gt2pt",
        help="call phenotypes from star alleles",
    )
    gt2pt_parser.add_argument(
        "gt",
        help="genotype file",
    )
    gt2pt_parser.add_argument(
        "-o",
        metavar="FILE",
        help="output to FILE [stdout]",
    )

    return parser

def output(fn, result):
    if fn:
        with open(fn, "w") as f:
            f.write(result)
    else:
        sys.stdout.write(result)

def main():
    start_time = timeit.default_timer()
    parser = get_parser()
    args = parser.parse_args()

    logger = logging.getLogger(__name__)
    logger.info(f"PyPGx v{__version__}")
    logger.info(f"""Command: '{" ".join(sys.argv)}'""")

    if args.tool == "pgkb":
        result = pgkb(args.t)
        output(args.o, result)

    elif args.tool == "report":
        result = report(args.gt)
        output(args.o, result)

    elif args.tool == "sdf2gdf":
        result = sdf2gdf(args.sdf, args.id)
        output(args.o, result)

    elif args.tool == "bam2sdf":
        result = bam2sdf(args.gb, args.tg, args.cg, args.bam)
        output(args.o, result)

    elif args.tool == "bam2gdf":
        bam = _get_bam_list(args.bam, args.bd, args.bl)
        result = bam2gdf(args.gb, args.tg, args.cg, bam)
        output(args.o, result)

    elif args.tool == "minivcf":
        result = minivcf(args.vcf, args.region)
        output(args.o, result)

    elif args.tool == "merge":
        result = merge(args.vcf, args.r)
        output(args.o, result)
    
    elif args.tool == "summary":
        result = summary(args.gt)
        output(args.o, result)

    elif args.tool == "meta":
        result = meta(args.sf)
        output(args.o, result)

    elif args.tool == "compare":
        result = compare(args.gt)
        output(args.o, result)

    elif args.tool == "remap":
        remap(args.conf)

    elif args.tool == "fq2bam":
        fq2bam(args.conf)

    elif args.tool == "sges":
        sges(args.conf)

    elif args.tool == "sgep":
        sgep(args.conf)

    elif args.tool == "sgea":
        sgea(args.conf)

    elif args.tool == "cpa":
        result = cpa(args.rdata)
        output(args.o, result)

    elif args.tool == "plotcov":
        plotcov(args.sdf, args.out)

    elif args.tool == "check":
        check(args.star, args.snp)

    elif args.tool == "liftover":
        result = liftover(args.star, args.snp, args.tg)
        output(args.o, result)

    elif args.tool == "peek":
        result = peek(args.vcf)
        output(args.o, result)

    elif args.tool == "snp":
        result = snp(args.vcf, args.pair)
        output(args.o, result)

    elif args.tool == "bam2vcf":
        bam = _get_bam_list(args.bam, args.bd, args.bl)
        result = bam2vcf(args.gb, args.tg, args.fa, bam)
        output(args.o, result)

    elif args.tool == "genotype":
        bam = _get_bam_list(args.bam, args.bd, args.bl)
        genotype(args.fa, args.dt, args.gb, args.tg, args.out, bam, args.cg)

    elif args.tool == "gt2pt":
        result = gt2pt(args.gt)
        output(args.o, result)

    else:
        pass

    stop_time = timeit.default_timer()
    elapsed_time = str(
        datetime.timedelta(seconds=(stop_time - start_time))
    ).split(".")[0]

    logger.info(f"Elapsed time: {elapsed_time}")
    logger.info("PyPGx finished")

if __name__ == "__main__":
    main()
