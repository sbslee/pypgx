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
    """
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
