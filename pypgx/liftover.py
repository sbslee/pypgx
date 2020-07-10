import os

from .common import get_gene_table

from .sglib import (
    read_snp_table,
    build_snpdb,
    read_star_table,
    build_stardb,
)

def liftover(star: str, snp: str, tg: str) -> str:
    """
    Convert variants in SNP table from hg19 to hg38.

    Returns:
        str: Result file.

    Args:
        star (str): Star allele table file.
        snp (str): SNP table file.
        tg (str): Target gene.
    """

    gene_table = get_gene_table()

    snp_table = read_snp_table(snp, gene_table)

    snpdb = build_snpdb(tg, "hg19", snp_table)

    star_table = read_star_table(star)

    stardb = build_stardb(tg, "hg19", star_table, snpdb)

    result = ""

    for k, v in stardb.items():
        core = []
        for snp in v.core:
            s = "{}:{}>{}".format(
                snp.data["hg38_pos"],
                snp.data["wt_allele"],
                snp.data["var_allele"]
            )
            core.append(s)
        if not core:
            core = ["."]
        tag = []
        for snp in v.tag:
            s = "{}:{}>{}".format(
                snp.data["hg38_pos"],
                snp.data["wt_allele"],
                snp.data["var_allele"]
            )
            tag.append(s)
        if not tag:
            tag = ["."]
        line = ",".join(core) + "\t" + ",".join(tag)
        result += line + "\n"

    return result