import os
import logging
from typing import Dict

import pysam

from .sglib import (
    read_gene_table,
    read_snp_table,
    build_snpdb,
    read_star_table,
    build_stardb,
    StarAllele,
)

LINE_BREAK1 = "-" * 70
LINE_BREAK2 = "*" * 70

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def sm_tag(bam: str) -> str:
    """
    Extract SM tag from BAM file.

    Returns:
        str: SM tag.

    Args:
        bam (str): BAM file.
    """

    header = pysam.view("-H", bam).strip().split("\n")

    l = []

    for line in header:
        fields = line.split("\t")
        if "@RG" == fields[0]:
            for field in fields:
                if "SM:" in field:
                    l.append(field.replace("SM:", ""))

    l = list(set(l))

    if not l:
        raise ValueError(f"SM tag not found: {bam}")

    if len(l) > 1:
        logger.warning(
            f"Multiple SM tags found (will return the first one): {bam}")
        result = l[0]
    else:
        result = l[0]

    return result

def is_chr(bam: str) -> bool:
    """
    Check whether SN tags in BAM file contain "chr" string.

    Returns:
        bool: True if found.

    Args:
        bam (str): BAM file.
    """

    header = pysam.view("-H", bam).strip().split("\n")

    l = []

    for line in header:
        fields = line.split("\t")
        if "@SQ" == fields[0]:
            for field in fields:
                if "SN:" in field:
                    l.append(field.replace("SN:", ""))

    return any(["chr" in x for x in l])

def get_stardb(tg: str, gb: str) -> Dict[str, StarAllele]:
    """
    Get StarAllele database.

    Returns:
        dict[str, StarAllele]: StarAllele objects.

    Args:
        tg (str): Target gene.
        gb (str): Genome build.

    Examples:

        >>> stardb = get_stardb("cyp2d6", "hg19")
        >>> print(stardb["*2"].score)
        1.0
    """

    p = os.path.dirname(__file__)
    gene_table = read_gene_table(f"{p}/resources/sg/gene_table.txt")
    snp_table = read_snp_table(f"{p}/resources/sg/snp_table.txt", gene_table)
    snpdb = build_snpdb(tg, gb, snp_table)
    star_table = read_star_table(f"{p}/resources/sg/star_table.txt")

    return build_stardb(tg, gb, star_table, snpdb)

def get_gene_table() -> Dict[str, Dict[str, str]]:
    """
    Get gene table object.

    Returns:
        dict[str, dict[str, str]]: Gene table object.
    """

    p = os.path.dirname(__file__)
    return read_gene_table(f"{p}/resources/sg/gene_table.txt")