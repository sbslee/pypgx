import os
import subprocess
from typing import Optional, List
from .common import get_target_region, is_chr, temp_env, bam_getter

def _run_haplotypecaller(
        fasta_file,
        input_file,
        gvcf_file,
        target_region
    ):
    command = [
        "gatk", "HaplotypeCaller",
        "-R", fasta_file,
        "--emit-ref-confidence", "GVCF",
        "-I", input_file,
        "-O", gvcf_file,
        "-L", target_region,
    ]

    subprocess.run(command, check=True)

def _run_combinegvcfs(
        fasta_file,
        target_region,
        gvcf_files,
        vcf_file
    ):
    command = [
        "gatk", "CombineGVCFs",
        "-R", fasta_file,
        "-L", target_region,
        "-O", vcf_file,
    ]

    for gvcf_file in gvcf_files:
        command += [
            "--variant", gvcf_file
        ]

    subprocess.run(command, check=True)

def _run_genotypegvcfs(
        fasta_file,
        dbsnp_file,
        target_region,
        vcf_file1,
        vcf_file2
    ):
    command = [
        "gatk", "GenotypeGVCFs",
        "-R", fasta_file,
        "-O", vcf_file2,
        "-L", target_region,
        "--variant", vcf_file1,
    ]

    if dbsnp_file:
        command += ["-D", dbsnp_file]

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
        temp_dir: Optional[str] = None,
        temp_path: str = None,
        input_files: List[str] = None,
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
            Target gene (e.g. 'cyp2d6') or region 
            (e.g.‘chr22:42512500-42551883’).
        output_file (str):
            Output VCF file.
        genome_build (str):
            Genome build ('hg19' or 'hg38').
        bam_file (list[str]):
            List of input BAM files.
        bam_dir (str, optional):
            Any BAM files in this directory will be used as input.
        bam_list (str, optional):
            List of BAM files, one file per line.
        dbsnp_file (str, optional):
            dbSNP VCF file used by GATK to add rs numbers.
        temp_dir (str, optional):
            Temporary files will be written to this directory (default: /tmp).
        temp_path (str):
            Automatically determined by @temp_env.
        input_files (list[str]):
            Automatically determined by @bam_getter.
    """
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
                target_region
            )

        _run_combinegvcfs(
            fasta_file,
            target_region,
            gvcf_files,
            f"{temp_path}/combined.g.vcf"
        )

        _run_genotypegvcfs(
            fasta_file,
            dbsnp_file,
            target_region,
            f"{temp_path}/combined.g.vcf",
            f"{temp_path}/combined.joint.vcf"
        )

        _run_variantfiltration(
            fasta_file,
            target_region,
            f"{temp_path}/combined.joint.vcf",
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
