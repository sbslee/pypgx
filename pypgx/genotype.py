from typing import Optional, List
from tempfile import TemporaryDirectory
from subprocess import run

from .bam2vcf import bam2vcf
from .bam2gdf import bam2gdf

def genotype(
        fa: str,
        sg: str,
        dt: str,
        gb: str,
        tg: str,
        out: str,
        bam: List[str],
        cg: Optional[str] = None,
    ) -> None:
    """Call star alleles from BAM file(s).

    This command runs the Stargazer genotyping pipeline without the need for 
    Sun Grid Engine (SGE). It uses the ``bam2vcf`` tool to create the input 
    VCF file (which is essentially a wrapper of, and therefore requires, the 
    BCFtools program) and the ``bam2gdf`` tool to create the input GDF file. 
    It then runs the Stargazer program to perform genotype analysis.

    In order to detect strctural variation, Stargazer needs read depth data 
    (i.e. a GDF file) for copy number analysis. Providing the optional 
    argument ``--cg`` will generate a GDF file.

    Args:
        fa (str):
            Reference FASTA file.
        sg (str):
            Stargazer program.
        dt (str):
            Sequencing data type. Use 'wgs' for whole 
            genome sequencing data and 'ts' for targeted 
            sequencing data.
        gb (str):
            Genome build ('hg19' or 'hg38').
        tg (str):
            Target gene (e.g. 'cyp2d6').
        out (str):
            Output project directory.
        bam (list[str]):
            List of BAM files.
        cg (str, optional):
            Control gene (e.g. 'vdr') or 
            region (e.g. ‘chr12:48232319-48301814’).

    .. note::

        BCFtools and Stargazer must be pre-installed.
    """
    tmp_dir = TemporaryDirectory()
    td = tmp_dir.name

    vcf = bam2vcf(gb, tg, fa, bam)

    with open(f"{td}/pypgx.vcf", "w") as f:
        f.write(vcf)

    if cg:
        gdf = bam2gdf(gb, tg, cg, bam)

        with open(f"{td}/pypgx.gdf", "w") as f:
            f.write(gdf)

    command = [
        "python3", sg,
        dt,
        gb,
        tg,
        f"{td}/pypgx.vcf",
        out,
    ]

    if cg:
        command += [
            "--cg", cg,
            "--gdf", f"{td}/pypgx.gdf",
        ]

    run(command)

    tmp_dir.cleanup()
