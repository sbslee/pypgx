from .phenotyper import phenotyper

def gt2pt(gt_file: str,
          **kwargs) -> str:
    """Convert a genotype file to phenotypes.

    This command is just a wrapper for the ``phenotyper`` module.

    Returns:
        str: Phenotype calls.

    Args:
        gt_file (str): Genotype file from Stargazer (``genotype.txt``).
    """

    result = []

    with open(gt_file) as f:
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