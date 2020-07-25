import subprocess
from typing import Optional, List
from .common import temp_env, bam_getter
from .bam2vcf import bam2vcf
from .bam2gdf import bam2gdf

@bam_getter
@temp_env
def bam2gt(
        snp_caller: str,
        fasta_file: str,
        target_gene: str,
        genome_build: str,
        data_type: str,
        proj_dir: str,
        bam_file: List[str],
        bam_dir: Optional[str] = None,
        bam_list: Optional[str] = None,
        control_gene: Optional[str] = None,
        dbsnp_file: Optional[str] = None,
        temp_dir: Optional[str] = None,
        temp_path: str = None,
        input_files: List[str] = None,
    ) -> None:
    """Convert BAM files to a genotype file.

    This command runs the entire genotyping pipeline for BAM files, 
    without the need for Sun Grid Engine (SGE). Under the hood, it 
    uses the ``bam2vcf`` command to create the input VCF file and 
    the ``bam2gdf`` command to create the input GDF file. It then 
    performs genotype analysis using the Stargazer program.

    In order to detect strctural variation, Stargazer requires read 
    depth data (i.e. a GDF file) for copy number analysis. Providing 
    the optional argument ``--control_gene`` will generate a GDF file. 
    If this argument is not provided, Stargazer will run as VCF-only mode.

    Args:
        snp_caller (str):
            SNP caller (‘gatk’ or ‘bcftools’).
        fasta_file (str):
            Reference FASTA file.
        target_gene (str):
            Target gene (e.g. ‘cyp2d6’) or 
            region (e.g.‘chr22:42512500-42551883’).
        genome_build (str):
            Genome build ('hg19' or 'hg38').
        data_type (str):
            Input data type ('wgs' or 'ts').
        proj_dir (str):
            Output project directory.
        bam_file (list[str]):
            List of input BAM files.
        bam_dir (str, optional):
            Any BAM files in this directory will be used as input.
        bam_list (str, optional):
            List of BAM files, one file per line.
        control_gene (str, optional):
        dbsnp_file (str, optional):
            dbSNP VCF file used by GATK to add rs numbers.
        temp_dir (str, optional):
            Temporary files will be written to this directory (default: /tmp).
        temp_path (str):
            Automatically determined by @temp_env.
        input_files (list[str]):
            Automatically determined by @bam_getter.
    """
    # Create the input VCF file.
    bam2vcf(
        snp_caller,
        fasta_file,
        target_gene,
        f"{temp_path}/pypgx.vcf",
        genome_build,
        bam_file = input_files
    )

    # Create the input GDF file.
    if control_gene:
        gdf = bam2gdf(
            genome_build,
            target_gene,
            control_gene,
            input_files
        )

        with open(f"{temp_path}/pypgx.gdf", "w") as f:
            f.write(gdf)

    # Run Stargazer.
    command = [
        "stargazer",
        data_type,
        genome_build,
        target_gene,
        f"{temp_path}/pypgx.vcf",
        proj_dir,
    ]

    if control_gene:
        command += [
            "--cg", control_gene,
            "--gdf", f"{temp_path}/pypgx.gdf",
        ]

    subprocess.run(command, check=True)
