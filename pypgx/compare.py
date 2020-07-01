import copy
from typing import List

def read_gt(x):
    genotype = []
    with open(x) as f:
        header = next(f).strip().split("\t")
        i_hap1_main = header.index("hap1_main")
        i_hap2_main = header.index("hap2_main")
        for line in f:
            fields = line.strip().split("\t")
            hap1_main = fields[i_hap1_main]
            hap2_main = fields[i_hap2_main]
            genotype.append([hap1_main, hap2_main])
    return genotype

def compare(gt: List[str]) -> str:
    """
    Compare genotype files.

    Returns:
        str: Result file.

    Args:
        gt (list[str]): Genotype file.
    """

    truth_file = gt[0]
    truth_genotype = read_gt(truth_file)

    n_total = len(truth_genotype)

    test_files = gt[1:]
    test_genotypes = []

    for test_file in test_files:
        test_genotypes.append(read_gt(test_file))

    data = []

    for test_genotype in test_genotypes:
        n_incorrect = 0
        for i in range(n_total):
            truth = copy.deepcopy(truth_genotype[i])
            test = test_genotype[i]

            for j in [0, 1]:
                if test[j] in truth:
                    truth.remove(test[j])

            n_incorrect += len(truth)

        n_correct = n_total * 2 - n_incorrect
        p_correct = n_correct / (n_total * 2)

        data.append((n_incorrect, n_correct, p_correct))

    temp = []
    temp.append(["name", "file", "disc", "cord", "perc"])
    temp.append(["truth", truth_file, 0, n_total * 2, 1.0])
    for i, x in enumerate(data):
        temp.append([f"test{i+1}", test_files[i], x[0], x[1], x[2]])

    result = ""

    for l in temp:
        result += "\t".join([
            "{0:.4f}".format(x) if isinstance(x, float) else str(x)
            for x in l]) + "\n"

    return result