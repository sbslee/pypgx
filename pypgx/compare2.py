from typing import Optional

def compare2(
        truth_file: str,
        test_file: str,
        sample_map: Optional[str] = None,
        **kwargs
    ) -> str:
    """Compare two genotype files.

    Returns:
        str: Genotype comparison.

    Args:
        truth_file (str):
            Truth genotype file.
        test_file (str):
            Test genotype file.
        sample_map (str, optional):
            Tab-delimited text file for sample name mapping.
    """

    mapping = {}

    with open(sample_map) as f:
        for line in f:
            fields = line.strip().split("\t")
            name1 = fields[0]
            name2 = fields[1]
            mapping[name1] = name2

    truth_data = {}

    with open(truth_file) as f:
        next(f)
        for line in f:
            fields = line.strip().split("\t")
            name = fields[1]
            hap1_main = fields[3]
            hap2_main = fields[4]
            genotype = hap1_main + "/" + hap2_main
            truth_data[name] = genotype

    test_data = {}

    with open(test_file) as f:
        next(f)
        for line in f:
            fields = line.strip().split("\t")
            name = fields[1]
            hap1_main = fields[3]
            hap2_main = fields[4]
            genotype = hap1_main + "/" + hap2_main
            test_data[name] = genotype

    result = ""

    for name1, genotype1 in truth_data.items():
        name2 = mapping[name1]
        genotype2 = test_data[name2]
        fields = [name1, genotype1, name2, genotype2, genotype1==genotype2]
        result += "\t".join([str(x) for x in fields]) + "\n"

    return result