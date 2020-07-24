import os
import subprocess
from tempfile import TemporaryDirectory

def _run_haplotypecaller(
        fasta_file,
        dbsnp_file,
        input_file,
        gvcf_file,
        target_region
    ):
    command = [
        "gatk", "HaplotypeCaller",
        "-R", fasta_file,
        "-D", dbsnp_file,
        "--emit-ref-confidence", "GVCF",
        "-I", input_file,
        "-O", gvcf_file,
        "-L", target_region,
    ]

    subprocess.run(command, check=True)

def _run_combinegvcfs(
        fasta_file,
        dbsnp_file,
        target_region,
        gvcf_files,
        vcf_file
    ):
    command = [
        "gatk", "CombineGVCFs",
        "-R", fasta_file,
        "-D", dbsnp_file,
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
        "-D", dbsnp_file,
        "-O", vcf_file2,
        "-L", target_region,
        "--variant", vcf_file1,
    ]

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

def bam2vcf2(
        snp_caller,
        fasta_file,
        dbsnp_file,
        target_region,
        output_file,
        bam_file,
        bam_dir = None,
        bam_list = None,
        temp_dir = None
    ):
    """Create a VCF file from BAM files.

    This command outputs a single- or multi-sample VCF file from one or 
    more input BAM files. The output VCF file will only contain variants
    within the target gene or region.
    """

    if temp_dir:
        temp_obj = None
        temp_path = os.path.realpath(temp_dir)
    else:
        temp_obj = TemporaryDirectory()
        temp_path = temp_obj.name

    if bam_file:
        input_files = bam_file
    elif bam_dir:
        bam_path = os.path.realpath(bam_dir)
        input_files = []
        for r, d, f in os.walk(bam_path):
            for x in f:
                if x.endswith("bam"):
                    input_files.append(f"{bam_path}/{x}")
    elif bam_list:
        input_files = []
        with open(bam_list) as f:
            for line in f:
                input_files.append(line.strip())
    else:
        input_files = None

    if not input_files:
        raise ValueError("No input BAM files found")

    if snp_caller == "gatk":

        gvcf_files = []

        for i, input_file in enumerate(input_files):
            gvcf_file = f"{temp_path}/{i}.g.vcf"
            gvcf_files.append(gvcf_file)

            _run_haplotypecaller(
                fasta_file,
                dbsnp_file,
                input_file,
                gvcf_file,
                target_region
            )

        _run_combinegvcfs(
            fasta_file,
            dbsnp_file,
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
        pass

    else:
        raise ValueError(f"Incorrect SNP caller: {snp_caller}")

    if temp_obj:
        temp_obj.cleanup()
