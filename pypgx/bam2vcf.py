import logging
from subprocess import Popen, PIPE, run
from typing import List
from tempfile import TemporaryDirectory

from .common import get_target_region, is_chr

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def bam2vcf(
        gb: str,
        tg: str,
        fasta: str,
        out: str,
        bam: List[str],
    ) -> None:
    """Create a VCF file from BAM file(s).

    This command outputs a single- or multi-sample VCF file from one or 
    more input BAM files. The output VCF file will only contain variants
    within the target gene or region. This is essentially a wrapper with
    certain parameters for various commands from the BCFtools program 
    (e.g. ``mpileup`` and ``call``). This means the called variants will be 
    already normalized and filtered, ready for the downstream genotype 
    analysis by the Stargazer program.

    Args:
        gb (str):
            Genome build ('hg19' or 'hg38').
        tg (str):
            Target gene (e.g. 'cyp2d6') or 
            region (e.g. 'chr22:42512500-42551883').
        fasta (str):
            Reference FASTA file.
        out (str):
            Output VCF file.
        bam (str):
            List of BAM files.

    .. note::

        BCFtools must be pre-installed.
    """
    tmp_dir = TemporaryDirectory()
    td = tmp_dir.name

    if ":" in tg:
        tr = tg.replace("chr", "")
    else:
        tr = get_target_region(tg, gb).replace("chr", "")

    t = [is_chr(x) for x in bam]

    if all(t):
        chr_str = "chr"
    elif not any(t):
        chr_str = ""
    else:
        raise ValueError("Mixed types of SN tags found.")

    command = [
        "bcftools", "mpileup",
        "-Ou",
         "-f", fasta,
         "-a", "AD",
         "-r", f"{chr_str}{tr}",
         "--max-depth", "1000",
    ] + bam

    p = Popen(command, stdout=PIPE)

    command = [
        "bcftools", "call",
        "-Oz",
        "-mv",
        "-o", f"{td}/calls.vcf.gz",
    ]

    run(command, stdin=p.stdout)

    p.stdout.close()

    run(["bcftools", "index", f"{td}/calls.vcf.gz"])

    run([
        "bcftools", "norm",
        f"{td}/calls.vcf.gz",
        "-Ob",
        "-f", fasta,
        "-o", f"{td}/calls.norm.bcf"
    ])

    run([
        "bcftools", "filter",
        f"{td}/calls.norm.bcf",
        "-Ov",
        "--IndelGap", "5",
        "-o", out,
    ])

    tmp_dir.cleanup()