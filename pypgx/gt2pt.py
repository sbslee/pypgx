from .phenotyper import phenotyper

def gt2pt(
        gt: str,
        **kwargs,
    ) -> str:
    """Call phenotypes from star alleles.

    This command is just a wrapper for the ``phenotyper`` module.

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

            result.append(phenotyper(gene, hap1, hap2))

    return "\n".join(result) + "\n"