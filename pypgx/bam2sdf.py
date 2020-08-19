import os
from typing import List

import pysam

from .common import logging, sm_tag, get_gene_table
from .sglib import sort_regions

logger = logging.getLogger(__name__)

def bam2sdf(
        genome_build: str,
        target_gene: str,
        control_gene: str,
        bam_file: List[str],
        **kwargs
    ) -> str:
    """
    Create SDF file from BAM file(s).

    Returns:
        str: SDF file.

    Args:
        genome_build (str): Genome build (hg19, hg38).
        target_gene (str): Target gene.
        control_gene (str): Control gene or region.
        bam_file (list[str]): BAM file(s).
    """

    gene_table = get_gene_table()

    targets = [k for k, v in gene_table.items() if v["type"] == "target"]

    if target_gene not in targets:
        raise ValueError(f"'{target_gene}' is not among target genes: {targets}")

    tr = gene_table[target_gene][f"{genome_build}_region"].replace("chr", "")

    if "chr" in control_gene or ":" in control_gene:
        cr = control_gene.replace("chr", "")

    else:
        controls = [k for k, v in gene_table.items() if v["control"] == "yes"]

        if control_gene not in controls:
            raise ValueError(f"'{control_gene}' is not among control genes: {controls}")

        cr = gene_table[control_gene][f"{genome_build}_region"].replace("chr", "")

    regions = sort_regions([tr, cr])

    # Get sample and sequence names from BAM headers.
    sm = []
    sn = []
    for x in bam_file:
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
        temp = pysam.depth("-a", "-Q", "1", "-r", f"{chr_str}{region}", *bam_file)
        result += temp

    return result
