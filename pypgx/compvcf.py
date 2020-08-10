from .sglib import VCFFile
from .common import get_logger

def read_sample_map(fn):
    result = []
    with open(fn) as f:
        for line in f:
            fields = line.strip().split("\t")
            result.append((fields[0], fields[1]))
    return result

def _get_summary(x):
    tn, tp, fn, fp = x
    con = (tn + tp) / (tn + tp + fn + fp)
    tpr = tp / (tp + fn) # sensitivity, recall, hit rate, or true positive rate (TPR)
    tnr = tn / (tn + fp) # specificity, selectivity or true negative rate (TNR)
    return ["{:.4f}".format(x) for x in [tpr, tnr, con]]

def compvcf(
        truth_file: str,
        test_file: str,
        sample_map: str,
        **kwargs
    ):
    """Compute the concordance between two VCF files.

    Returns:
        str: Genotype concordance.

    Args:
        truth_file (str):
            Truth VCF file.
        test_file (str):
            Test VCF file.
        sample_map (str, optional):
            Tab-delimited text file with two columns representing the truth 
            and test sample names.
    """
    logger = get_logger()

    mapping = read_sample_map(sample_map)

    vcf1 = VCFFile(truth_file)
    vcf1.read(tidy=True)
    logger.info(f"Total number of SNPs in truth VCF: {len(vcf1.data)}")
    vcf1.multiallelic_filter()
    logger.info(f"Total number of SNPs in truth VCF after multiallelic filter: {len(vcf1.data)}")
    vcf1.missing_filter()
    logger.info(f"Total number of SNPs in truth VCF after missing filter: {len(vcf1.data)}")

    vcf2 = VCFFile(test_file)
    vcf2.read(tidy=True)
    logger.info(f"Total number of SNPs in test VCF: {len(vcf2.data)}")
    vcf2.multiallelic_filter()
    logger.info(f"Total number of SNPs in test VCF after multiallelic filter: {len(vcf2.data)}")
    vcf2.missing_filter()
    logger.info(f"Total number of SNPs in truth VCF after missing filter: {len(vcf2.data)}")

    snp1 = ["{}:{}:{}:{}".format(x[0], x[1], x[3], x[4]) for x in vcf1.data]
    snp2 = ["{}:{}:{}:{}".format(x[0], x[1], x[3], x[4]) for x in vcf2.data]

    overlapping = sorted(list(set(snp1) & set(snp2)))

    logger.info(f"Total number of overlapping SNPs: {len(overlapping)}")

    for i in reversed(range(len(vcf1.data))):
        v = vcf1.data[i]
        x = "{}:{}:{}:{}".format(v[0], v[1], v[3], v[4])
        if x not in overlapping:
            del vcf1.data[i]

    for i in reversed(range(len(vcf2.data))):
        v = vcf2.data[i]
        x = "{}:{}:{}:{}".format(v[0], v[1], v[3], v[4])
        if x not in overlapping:
            del vcf2.data[i]

    dat = []

    for i in range(len(mapping)):
        name1 = mapping[i][0]
        name2 = mapping[i][1]

        x1 = vcf1.header.index(name1)
        x2 = vcf2.header.index(name2)

        # tn, tp, fn, fp
        counts = {
            "snv": [0, 0, 0, 0],
            "indel": [0, 0, 0, 0]
        }

        for j in range(len(vcf1.data)):
            v1 = vcf1.data[j]
            v2 = vcf2.data[j]

            if len(v1[3]) == 1 and len(v1[4]) == 1:
                type = "snv"
            else:
                type = "indel"

            ac1 = sum([int(x) for x in v1[x1].split(":")[0].replace("|", "/").split("/")])
            ac2 = sum([int(x) for x in v2[x2].split(":")[0].replace("|", "/").split("/")])

            if ac1 == ac2 == 0:
                counts[type][0] += 1
            elif ac1 == ac2:
                counts[type][1] += 1
            elif ac1 > ac2:
                counts[type][2] += 1
            elif ac1 < ac2:
                counts[type][3] += 1

        snv = counts["snv"]
        indel = counts["indel"]
        all = [x + y for x, y in zip(counts["snv"], counts["indel"])]

        row = [name1, name2] + snv + _get_summary(snv) + indel + _get_summary(indel) + all + _get_summary(all)

        dat.append(row)

    result = ""

    headers = [
        "name1", "name2", "snv_tn", "snv_tp", "snv_fn", "snv_fp", "snv_tpr", "snv_tnr", "snv_con",
        "indel_tn", "indel_tp", "indel_fn", "indel_fp", "indel_tpr", "indel_tnr", "indel_con",
        "all_tn", "all_tp", "all_fn", "all_fp", "all_tpr", "all_tnr", "all_con",
    ]

    result += "\t".join(headers) + "\n"

    for row in dat:
        result += "\t".join([str(x) for x in row]) + "\n"

    return result
