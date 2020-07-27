import os
import subprocess
from typing import Optional, List
from .common import get_target_region, is_chr, temp_env, bam_getter

def _run_haplotypecaller(
        fasta_file,
        input_file,
        gvcf_file,
        target_region,
        java_options
    ):
    command = [
        "gatk", "HaplotypeCaller",
        "-R", fasta_file,
        "--emit-ref-confidence", "GVCF",
        "-I", input_file,
        "-O", gvcf_file,
        "-L", target_region,
        "--QUIET",
    ]

    if java_options:
        command += ["--java-options", java_options]

    subprocess.run(command, check=True)

def _run_genomicsdbimport(
        target_region,
        gvcf_files,
        datastore
    ):
    command = [
        "gatk", "GenomicsDBImport",
        "--intervals", target_region,
        "--genomicsdb-workspace-path", datastore,
        "--QUIET",
    ]

    for gvcf_file in gvcf_files:
        command += [
            "-V", gvcf_file
        ]

    subprocess.run(command, check=True)

def _run_genotypegvcfs(
        fasta_file,
        dbsnp_file,
        datastore,
        vcf_file,
        java_options
    ):
    command = [
        "gatk", "GenotypeGVCFs",
        "-R", fasta_file,
        "-V", datastore,
        "-O", vcf_file,
        "--QUIET",
    ]

    if dbsnp_file:
        command += ["-D", dbsnp_file]

    if java_options:
        command += ["--java-options", java_options]

    subprocess.run(command, check=True)

def _run_variantfiltration(
        fasta_file,
        target_region,
        vcf_file,
        output_file
    ):
    command = [
        "gatk", "VariantFiltration",
        "-R", fasta_file,
        "-L", target_region,
        "-O", output_file,
        "--filter-expression", "QUAL <= 50.0",
        "--filter-name", "QUALFilter",
        "--variant", vcf_file,
        "--QUIET",
    ]

    subprocess.run(command, check=True)

def _run_mpileup(
        fasta_file,
        target_region,
        input_files,
        vcf_file
    ):
    command = [
        "bcftools", "mpileup",
        "-Ou",
        "-f", fasta_file,
        "-a", "AD",
        "-r", target_region,
        "--max-depth", "1000",
        "-o", vcf_file,
    ] + input_files

    subprocess.run(command, check=True)

def _run_call(
        vcf_file1,
        vcf_file2
    ):
    command = [
        "bcftools", "call",
        vcf_file1,
        "-Oz",
        "-mv",
        "-o", vcf_file2
    ]

    subprocess.run(command, check=True)

def _run_index(vcf_file):
    command = [
        "bcftools", "index",
        vcf_file,
    ]

    subprocess.run(command, check=True)

def _run_norm(
        fasta_file,
        vcf_file1,
        vcf_file2
    ):
    command = [
        "bcftools", "norm",
        vcf_file1,
        "-Ob",
        "-f", fasta_file,
        "-o", vcf_file2,
    ]

    subprocess.run(command, check=True)

def _run_filter(
        vcf_file1,
        output_file
    ):
    command = [
        "bcftools", "filter",
        vcf_file1,
        "-Ov",
        "--IndelGap", "5",
        "-o", output_file,
    ]

    subprocess.run(command, check=True)

@bam_getter
@temp_env
def bam2vcf(
        snp_caller: str,
        fasta_file: str,
        target_gene: str,
        output_file: str,
        genome_build: str,
        bam_file: List[str],
        bam_dir: Optional[str] = None,
        bam_list: Optional[str] = None,
        dbsnp_file: Optional[str] = None,
        java_options: Optional[str] = None,
        temp_dir: Optional[str] = None,
        **kwargs
    ) -> None:
    """Create a VCF file from BAM files.

    This command creates a single- or multi-sample VCF file from one or 
    more input BAM files. The output VCF file will only contain variants 
    within the target gene or region. The command is essentially a wrapper 
    for the Genome Analysis Toolkit (GATK) and the BCFtools program with 
    pre-specified parameters. This means the called variants will be 
    already normalized and filtered, ready for the downstream genotype 
    analysis by the Stargazer program.

    Args:
        snp_caller (str):
            SNP caller ('gatk' or 'bcftools').
        fasta_file (str):
            Reference FASTA file.
        target_gene (str):
            Target gene (e.g. 'cyp2d6') or 
            region (e.g.‘chr22:42512500-42551883’).
        output_file (str):
            Write output to this file.
        genome_build (str):
            Genome build ('hg19' or 'hg38').
        bam_file (list[str]):
            Input BAM files.
        bam_dir (str, optional):
            Use all BAM files in this directory as input.
        bam_list (str, optional):
            List of input BAM files, one file per line.
        dbsnp_file (str, optional):
            dbSNP VCF file used by GATK to add rs numbers.
        java_options (str, optional):
            Java-specific arguments for GATK (e.g. '-Xmx4G').
        temp_dir (str, optional):
            Temporary files will be written to this directory.

    .. warning::
        GATK and/or BCFtools must be pre-installed.

    .. note::
        Generally, GATK is more accurate but much slower than BCFtools. 
        For instance, SNP calling for 70 WGS samples for the CYP2D6 gene 
        takes 19 min with the ``gatk`` caller but only 2 min with the 
        ``bcftools`` caller. Therefore, if you have many samples and you do 
        not have access to Sun Grid Engine (SGE) for parallelism, we 
        recommend that you use ``bcftools``. If you have access to SGE and 
        you want to use GATK instead of BCFtools, please check other 
        SGE-based commands in PyPGx (e.g. ``sgep``).
    """
    # Parse keyward arguments from the decorators.
    temp_path = kwargs["temp_path"]
    input_files = kwargs["input_files"]

    # Pick the chromosome string.
    _ = [is_chr(x) for x in input_files]

    if all(_):
        chr_str = "chr"
    elif not any(_):
        chr_str = ""
    else:
        raise ValueError("Mixed types of SN tags found.")    

    # Get the study region.
    if ":" in target_gene:
        target_region = chr_str + target_gene.replace("chr", "")
    else:
        target_region = chr_str + get_target_region(
            target_gene, genome_build).replace("chr", "")

    # Run the selected SNP caller.
    if snp_caller == "gatk":

        gvcf_files = []

        for i, input_file in enumerate(input_files):
            gvcf_file = f"{temp_path}/{i}.g.vcf"
            gvcf_files.append(gvcf_file)

            _run_haplotypecaller(
                fasta_file,
                input_file,
                gvcf_file,
                target_region,
                java_options
            )

        _run_genomicsdbimport(
            target_region,
            gvcf_files,
            f"{temp_path}/datastore"
        )

        _run_genotypegvcfs(
            fasta_file,
            dbsnp_file,
            f"gendb://{temp_path}/datastore",
            f"{temp_path}/.joint.vcf",
            java_options
        )

        _run_variantfiltration(
            fasta_file,
            target_region,
            f"{temp_path}/.joint.vcf",
            output_file
        )

    elif snp_caller == "bcftools":
        _run_mpileup(
            fasta_file,
            target_region,
            input_files,
            f"{temp_path}/uncompressed.bcf"
        )

        _run_call(
            f"{temp_path}/uncompressed.bcf",
            f"{temp_path}/calls.vcf.gz"
        )

        _run_index(
            f"{temp_path}/calls.vcf.gz"
        )

        _run_norm(
            fasta_file,
            f"{temp_path}/calls.vcf.gz",
            f"{temp_path}/calls.norm.bcf"
        )

        _run_filter(
            f"{temp_path}/calls.norm.bcf",
            output_file
        )

    else:
        raise ValueError(f"Incorrect SNP caller: {snp_caller}")
