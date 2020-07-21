from .common import get_stardb

def _phenotype_cyp2c19(stardb, hap1, hap2):
    as1 = _hap2as(stardb, hap1)
    as2 = _hap2as(stardb, hap2)
    total = as1 + as2
    if total < 0:
        result = "unknown_metabolizer"
    elif total == 0:
        result = "poor_metabolizer"
    elif 0 < total <= 1.25:
        result = "intermediate_metabolizer"
    elif 1.25 < total <= 2:
        result = "normal_metabolizer"
    elif 2 < total < 2.5:
        result = "rapid_metabolizer"
    elif total >= 2.5:
        result = "ultrarapid_metabolizer"
    else:
        result = "undetermined_metabolizer"
    return result

def _phenotype_cyp2d6(stardb, hap1, hap2):
    as1 = _hap2as(stardb, hap1)
    as2 = _hap2as(stardb, hap2)
    total = as1 + as2
    if total < 0:
        result = "unknown_metabolizer"
    elif total == 0:
        result = "poor_metabolizer"
    elif 0 < total <= 1:
        result = "intermediate_metabolizer"
    elif 1 < total <= 2.25:
        result = "normal_metabolizer"
    elif total > 2.25:
        result = "ultrarapid_metabolizer"
    else:
        result = "undetermined_metabolizer"
    return result

def _hap2as(stardb, hap):
    result = 0
    for sa in hap.split("+"):
        if "x" in sa:
            n = int(sa.split("x")[1])
            name = sa.split("x")[0]
            result += stardb[name].score * n
        else:
            result += stardb[sa].score
    return result

def phenotype(
        gene: str,
        hap1: str,
        hap2: str
    ) -> str:
    """Predict a phenotype from two haplotype calls.

    Returns:
        str: Phenotype call.

    Args:
        gene (str): Target gene.
        hap1 (str): 1st haplotype call.
        hap2 (str): 2nd haplotype call.
    """
    stardb = get_stardb(gene, "hg19")

    if gene == "cyp2c19":
        result = _phenotype_cyp2c19(stardb, hap1, hap2)

    elif gene == "cyp2d6":
        result = _phenotype_cyp2d6(stardb, hap1, hap2)

    else:
        result = f"There are no known phenotypes for the gene: {gene}"

    return result

def gt2pt(gt: str) -> str:
    """Predict phenotypes from star allele calls.

    This method is just a wrapper for reading the genotype file. 
    The actual phenotyping is performed by the ``phenotype`` method.

    Returns:
        str: Phenotype calls.

    Args:
        gt (str): Genotype file from Stargazer.
    """
    result = []

    with open(gt) as f:
        header = next(f).strip().split("\t")

        i1 = header.index("gene")
        i2 = header.index("hap1_main")
        i3 = header.index("hap2_main")

        for line in f:
            fields = line.strip().split("\t")
            gene = fields[i1]
            hap1 = fields[i2]
            hap2 = fields[i3]

            result.append(phenotype(gene, hap1, hap2))

    return "\n".join(result) + "\n"