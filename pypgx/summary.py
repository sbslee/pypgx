import os

from .sglib import sort_star_names, read_phenotype_table

def summary(tg: str, gt: str) -> str:
    """
    Create summary file using Stargazer data.

    Returns:
        str: Summary file.

    Args:
        tg (str): Target gene.
        gt (str): Genotype file.
    """

    phenotype_table = read_phenotype_table(
        f"{os.path.dirname(__file__)}/resources/sg/phenotype_table.txt")

    star_dict = {}   # list all observed star alleles (haplotypes)
    asc_dict = {}    # how many observed alleles had AS=0, 0.5, ...
    hsc_dict = {}    # how many sample haplotypes had AS=0, 0.5, ...
    dsc_dict = {}    # how many samples had AS=0, 0.5, ...
    pt_dict = {}     # list all observed phenotypes
    sampsv_dict = {} # how many samples had zero SVs, one SV, or two SVs 
    sv_dict = {}     # list all observed SV calls

    s_sv = 0
    a_total = 0
    s_total = 0

    with open(gt) as f:
        next(f)
        for line in f:
            fields = line.strip().split("\t")
            hap1_main = fields[2]
            hap2_main = fields[3]
            hap1_score = fields[6]
            hap2_score = fields[7]
            dip_score = fields[8]
            phenotype = fields[9]
            dip_sv = fields[10]
            hap1_sv = fields[11]
            hap2_sv = fields[12]

            if hap1_main not in star_dict:
                star_dict[hap1_main] = [hap1_sv, hap1_score, 0]
            if hap2_main not in star_dict:
                star_dict[hap2_main] = [hap2_sv, hap2_score, 0]
            if dip_score not in dsc_dict:
                dsc_dict[dip_score]= 0
            if phenotype not in pt_dict:
                pt_dict[phenotype] = 0
            if dip_sv != "no_sv,no_sv":
                s_sv += 1

            if "unknown" in dip_sv:
                raise ValueError()
            else:
                sv_count = 2 - dip_sv.split(",").count("no_sv")

            if sv_count not in sampsv_dict:
                sampsv_dict[sv_count] = 0
            if hap1_score not in hsc_dict:
                hsc_dict[hap1_score] = 0
            if hap2_score not in hsc_dict:
                hsc_dict[hap2_score] = 0
            if hap1_sv not in sv_dict:
                sv_dict[hap1_sv] = 0
            if hap2_sv not in sv_dict:
                sv_dict[hap2_sv] = 0

            hsc_dict[hap1_score] += 1
            hsc_dict[hap2_score] += 1
            star_dict[hap1_main][2] += 1
            star_dict[hap2_main][2] += 1
            dsc_dict[dip_score] += 1
            pt_dict[phenotype] += 1
            sampsv_dict[sv_count] += 1
            sv_dict[hap1_sv] += 1
            sv_dict[hap2_sv] += 1

    a_total = sum([x[2] for x in star_dict.values()])
    s_total = int(a_total / 2)
    
    # Count allele scores
    for allele in star_dict:
        score = star_dict[allele][1]
        if score not in asc_dict:
            asc_dict[score] = 0
        asc_dict[score] += 1

    # Combine everything into a list of lists.
    temp = []
    temp.append(["type", "name", "sv", "sc", "n", "p"])
    temp.append(["samp", "total", ".", ".", s_total, s_total/s_total])
    temp.append(["samp", "sv", ".", ".", s_sv, s_sv/s_total])

    for sv_count in sorted(sampsv_dict):
        n = sampsv_dict[sv_count]
        p = n / s_total
        temp.append(["sampsv", str(sv_count), ".", ".", n, p])

    temp.append(["haps", "total", ".", ".", a_total, a_total/a_total])
    temp.append(["haps", "unique", ".", ".", len(star_dict), "."])

    for sv_type in sorted(sv_dict):
        n = sv_dict[sv_type]
        p = n / a_total
        temp.append(["sv", sv_type, ".", ".", n, p])

    for allele in sort_star_names(list(star_dict)):
        sv = star_dict[allele][0]
        s = star_dict[allele][1]
        n = star_dict[allele][2]
        p = n / a_total
        temp.append(["star", allele, sv, s, n, p])

    for allele_score in sorted(asc_dict):
        n = asc_dict[allele_score]
        p = n / len(star_dict)
        temp.append(["asc", allele_score, ".", ".", n, p])

    for hap_score in sorted(hsc_dict):
        n = hsc_dict[hap_score]
        p = n / a_total
        temp.append(["hsc", hap_score, ".", ".", n, p])

    for dip_score in sorted(dsc_dict):
        n = dsc_dict[dip_score]
        p = n / s_total
        temp.append(["dsc", dip_score, ".", ".", n, p])

    for x in sorted(pt_dict, key = list(phenotype_table[tg]).index):
        n = pt_dict[x]
        p = pt_dict[x] / s_total
        temp.append(["pt", x, ".", ".", n, p])

    result = ""

    for l in temp:
        result += "\t".join([
            "{0:.4f}".format(x) if isinstance(x, float) else str(x)
            for x in l]) + "\n"

    return result