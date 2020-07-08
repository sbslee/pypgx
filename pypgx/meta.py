import os
from typing import List

from .sglib import sort_star_names, read_phenotype_table

def _append1(d, i, name, count, p, sv, hap_score):
    if name not in d:
        d[name] = [sv, hap_score]
        for j in range(i):
            d[name].append(".")
            d[name].append(".")
    d[name].append(count)
    d[name].append(p)

def _append2(d, i):
    for name in d:
        if len(d[name]) == 2 * i + 2:
            d[name].append(".")
            d[name].append(".")

def meta(tg: str, sf: List[str]) -> str:
    """
    Create meta file from summary files.

    Returns:
        str: Meta file.

    Args:
        tg (str): Target gene.
        sf (list[str]): Summary file.
    """

    phenotype_table = read_phenotype_table(
        f"{os.path.dirname(__file__)}/resources/sg/phenotype_table.txt")

    dicts = {}
    header1 = ["type", "name", "sv", "hap_score"]

    for i in range(len(sf)):
        summary = sf[i]
        prefix = os.path.basename(summary)
        header1.append("n_" + prefix)
        header1.append("p_" + prefix)
        with open(summary) as f:
            header2 = next(f).strip().split("\t")
            for line in f:
                fields = line.strip().split("\t")
                type = fields[0]
                name = fields[1]
                sv = fields[2]
                hap_score = fields[3]
                count = fields[4]
                p = fields[5]
                if type not in dicts:
                    dicts[type] = {}
                _append1(dicts[type], i, name, count, p, sv, hap_score)
        for type in dicts:
            _append2(dicts[type], i)

    temp = []
    temp.append(header1)

    for type in dicts:
        if type == "samp":
            for x in ["total", "sv"]:
                temp.append(["samp", x] + dicts[type][x])

        elif type == "haps":
            for x in ["total", "unique"]:
                temp.append(["haps", x] + dicts[type][x])

        elif type == "star":
            for x in sort_star_names(list(dicts[type])):
                temp.append(["star", x] + dicts[type][x])

        elif type == "pt":
            for x in sorted(
                    list(dicts[type]),
                    key = list(phenotype_table[tg]).index
                ):
                temp.append(["pt", x] + dicts[type][x])

        else:
            for x in sorted(dicts[type]):
                temp.append([type, x] + dicts[type][x])

    result = ""

    for l in temp:
        result += "\t".join(l) + "\n"

    return result
