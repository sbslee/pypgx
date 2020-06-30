import os
from typing import List

from pypgx.common import sort_star_names, read_pt_table

def _append1(d, i, name, count, percentage, sv, hap_score):
    if name not in d:
        d[name] = [sv, hap_score]
        for j in range(i):
            d[name].append(".")
            d[name].append(".")
    d[name].append(count)
    d[name].append(percentage)

def _append2(d, i):
    for name in d:
        if len(d[name]) == 2 * i + 2:
            d[name].append(".")
            d[name].append(".")

def meta(tg: str, sum: List[str]) -> str:
    """
    Create meta file from summary files.

    Returns:
        str: Meta file.

    Args:
        tg (str): Target gene.
        sum (list[str]): Summary file.
    """

    dicts = {}
    header1 = ["type", "name", "sv", "hap_score"]

    for i in range(len(sum)):
        summary = sum[i]
        prefix = os.path.basename(summary)
        header1.append("N_" + prefix)
        header1.append("P_" + prefix)
        with open(summary) as f:
            header2 = next(f).strip().split("\t")
            for line in f:
                fields = line.strip().split("\t")
                type = fields[0]
                name = fields[1]
                sv = fields[2]
                hap_score = fields[3]
                count = fields[4]
                percentage = fields[5]
                if type not in dicts:
                    dicts[type] = {}
                _append1(dicts[type], i, name, count, percentage, sv, hap_score)
        for type in dicts:
            _append2(dicts[type], i)


    temp = []
    temp.append(header1)

    for type in dicts:
        if type == "samp":
            for subtype in ["total", "sv"]:
                temp.append(["samp", subtype] + dicts[type][subtype])

        elif type == "haps":
            for subtype in ["total", "unique"]:
                temp.append(["haps", subtype] + dicts[type][subtype])

        elif type == "star":
            for allele in sort_star_names(list(dicts[type])):
                temp.append(["star", allele] + dicts[type][allele])

        elif type == "pt":
            for phenotype in sorted(list(dicts[type]), key = list(read_pt_table()[tg]).index):
                temp.append(["pt", phenotype] + dicts[type][phenotype])

        else:
            for subtype in sorted(dicts[type]):
                temp.append([type, subtype] + dicts[type][subtype])

    result = ""

    for l in temp:
        result += "\t".join(l) + "\n"

    return result
