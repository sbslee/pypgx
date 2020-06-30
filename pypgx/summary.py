from pypgx.common import read_pt_table, sort_star_names

def summary(gt: str, tg: str) -> str:
    """
    Create summary file using data from Stargazer.

    Returns:
        str: Summary file.

    Args:
        gt (str): Stargazer genotype file.
        tg (str): Target gene.
    """

    pt_table = read_pt_table()
    
    alleles = {}
    allele_scores = {}
    hap_scores = {}
    dip_scores = {}
    phenotypes = {}
    sv_counts = {}
    sv_types = {}
    s_sv = 0

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

            if hap1_main not in alleles:
                alleles[hap1_main] = [hap1_sv, hap1_score, 0]
            if hap2_main not in alleles:
                alleles[hap2_main] = [hap2_sv, hap2_score, 0]
            if dip_score not in dip_scores:
                dip_scores[dip_score]= 0
            if phenotype not in phenotypes:
                phenotypes[phenotype] = 0
            if dip_sv != "no_sv,no_sv":
                s_sv += 1
            if "unknown" in dip_sv:
                raise ValueError()
            else:
                sv_count = 2 - dip_sv.split(",").count("no_sv")
            if sv_count not in sv_counts:
                sv_counts[sv_count] = 0
            if hap1_score not in hap_scores:
                hap_scores[hap1_score] = 0
            if hap2_score not in hap_scores:
                hap_scores[hap2_score] = 0
            if hap1_sv not in sv_types:
                sv_types[hap1_sv] = 0
            if hap2_sv not in sv_types:
                sv_types[hap2_sv] = 0
            hap_scores[hap1_score] += 1
            hap_scores[hap2_score] += 1
            alleles[hap1_main][2] += 1
            alleles[hap2_main][2] += 1
            dip_scores[dip_score] += 1
            phenotypes[phenotype] += 1
            sv_counts[sv_count] += 1
            sv_types[hap1_sv] += 1
            sv_types[hap2_sv] += 1
    a_total = sum([x[2] for x in alleles.values()])    
    s_total = int(a_total / 2)
    
    # Count allele scores
    for allele in alleles:
        score = alleles[allele][1]
        if score not in allele_scores:
            allele_scores[score] = 0
        allele_scores[score] += 1
        
    
    temp = []
    temp.append(["type", "name", "sv", "sc", "n", "perc"])
    temp.append(["samp", "total", ".", ".", s_total, s_total/s_total])
    temp.append(["samp", "sv", ".", ".", s_sv, s_sv/s_total])
    temp.append(["sampsv", str(sv_count), ".", ".", sv_counts[sv_count], sv_counts[sv_count] / s_total])
    temp.append(["haps", "total", ".", ".", a_total, a_total/a_total])
    temp.append(["haps", "unique", ".", ".", len(alleles), "."])
    for sv_type in sorted(sv_types):
        temp.append(["sv", sv_type, ".", ".", sv_types[sv_type], sv_types[sv_type] / a_total])

    for allele in sort_star_names(list(alleles)):
        sv = alleles[allele][0]
        hap_score = alleles[allele][1]
        n = alleles[allele][2]
        p = n / a_total
        temp.append(["star", allele, sv, hap_score, n, p])

    for allele_score in sorted(allele_scores):
        n = allele_scores[allele_score]
        p = n / len(alleles)
        temp.append(["asc", allele_score, ".", ".", n, p])

    for hap_score in sorted(hap_scores):
        n = hap_scores[hap_score]
        p = n / a_total
        temp.append(["hsc", hap_score, ".", ".", n, p])

    for dip_score in sorted(dip_scores):
        n = dip_scores[dip_score]
        p = n / s_total
        temp.append(["dsc", dip_score, ".", ".", n, p])

    for phenotype in sorted(phenotypes, key = list(pt_table[tg]).index):
        n = phenotypes[phenotype]
        p = phenotypes[phenotype] / s_total
        temp.append(["pt", phenotype, ".", ".", n, p])

    result = ""

    for l in temp:
        result += "\t".join([
            "{0:.4f}".format(x) if isinstance(x, float) else str(x) for x in l
        ]) + "\n"

    return result