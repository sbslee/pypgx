import logging
from subprocess import Popen, PIPE
from typing import List

from .common import get_target_region, is_chr

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def bam2vcf(
        gb: str,
        tg: str,
        fasta: str,
        bam: List[str],
    ) -> str:
    """Create a single- or multi-sample VCF file from one or more BAM files.

    Returns:
        str: A VCF file.

    Args:
        gb (str):
            Genome build (hg19 or hg38).
        tg (str):
            Target gene (e.g. 'cyp2d6') or 
            region (e.g. 'chr22:42512500-42551883').
        fasta (str):
            Reference FASTA file.
        bam (str):
            List of BAM files.

    .. note::

        This is essentially a wrapper for the ``mpileup`` and ``call`` 
        commands from the BCFtools program; therefore, you must have 
        BCFtools installed before running this method.
    """
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
        "bcftools", "mpileup", "-Ou",
         "-f", fasta,
         "-a", "AD",
         "-r", f"{chr_str}{tr}",
    ]

    logger.info("Command for bcftools-mpileup:")
    logger.info("    '{}'".format(" ".join(command + bam)))

    p1 = Popen(
        command + bam,
        stdout=PIPE,
        stderr=PIPE
    )

    command = ["bcftools", "call", "-mv", "-Ov"]

    logger.info("Command for bcftools-call:")
    logger.info("    '{} FILE'".format(" ".join(command)))

    p2 = Popen(
        command,
        stdin=p1.stdout, # p1.stdout is a real file
        stdout=PIPE,
        stderr=PIPE
    )

    e1 = p1.stderr.read()

    p1.stdout.close()
    p1.stderr.close()

    out, e2 = p2.communicate()

    logger.info("Message from bcftools-mpileup:")
    logger.info(f"    {e1}")

    logger.info("Message from bcftools-call:")
    logger.info(f"    {e2}")
 
    result = out.decode("utf-8")

    return result
