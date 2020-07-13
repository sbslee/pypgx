import os
from typing import List

import pysam

from .common import logging, sm_tag, get_gene_table
from .sglib import sort_regions

logger = logging.getLogger(__name__)

def bam2sdf(
        gb:str,
        tg: str,
        cg: str,
        bam: List[str]
    ) -> str:
    """
    Create SDF file from BAM file(s).

    Returns:
        str: SDF file.

    Args:
        gb (str): Genome build (hg19, hg38).
        tg (str): Target gene.
        cg (str): Control gene or region.
        bam (list[str]): BAM file(s).
    """

    gene_table = get_gene_table()

    targets = [k for k, v in gene_table.items() if v["type"] == "target"]

    if tg not in targets:
        raise ValueError(f"'{tg}' is not among target genes: {targets}")

    tr = gene_table[tg][f"{gb}_region"].replace("chr", "")

    if "chr" in cg or ":" in cg:
        cr = cg.replace("chr", "")

    else:
        controls = [k for k, v in gene_table.items() if v["control"] == "yes"]

        if cg not in controls:
            raise ValueError(f"'{cg}' is not among control genes: {controls}")

        cr = gene_table[cg][f"{gb}_region"].replace("chr", "")

    regions = sort_regions([tr, cr])

    # Get sample and sequence names from BAM headers.
    sm = []
    sn = []
    for x in bam:
        sm.append(sm_tag(x))

        result = pysam.view("-H", x).strip().split("\n")
        for line in result:
            fields = line.split("\t")
            if "@SQ" == fields[0]:
                for field in fields:
                    if "SN:" in field:
                        y = field.replace("SN:", "")
                        if y not in sn:
                            sn.append(y)

    logger.info(f"Sample IDs: {sm}")
    logger.info(f"Contigs: {sn}")

    # Determine whether the "chr" string should be used.
    if any(["chr" in x for x in sn]):
        chr_str = "chr"
    else:
        chr_str = ""

    result = ""

    for region in regions:
        temp = pysam.depth("-a", "-Q", "1", "-r", f"{chr_str}{region}", *bam)
        result += temp

    return result
