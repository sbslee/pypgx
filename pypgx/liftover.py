import copy

from .common import build_stardb

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

    star_db = build_stardb(tg, star, snp)

    result = ""

    for k, v in star_db.items():
        core = []
        for snp in v.core:
            s = "{}:{}>{}".format(snp.data["hg38_pos"], snp.data["wt_allele"], snp.data["var_allele"])
            core.append(s)
        if not core:
            core = ["."]
        tag = []
        for snp in v.tag:
            s = "{}:{}>{}".format(snp.data["hg38_pos"], snp.data["wt_allele"], snp.data["var_allele"])
            tag.append(s)
        if not tag:
            tag = ["."]
        line = ",".join(core) + "\t" + ",".join(tag)
        result += line + "\n"

    return result