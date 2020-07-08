import os

from .common import get_stardb
from .sglib import VCFFile, parse_vcf_fields, vcf2samples

def peek(vcf) -> str:
    """
    Find all possible star alleles from VCF file.

    Returns:
        str: Result file.

    Args:
        vcf (str): VCF file.
    """

    # Remove sample data from the VCF file.
    finalized_vcf = VCFFile(vcf)
    finalized_vcf.read()
    finalized_vcf.close()
    finalized_vcf.header = finalized_vcf.header[:9]

    stardb = get_stardb(finalized_vcf.tg, finalized_vcf.gb)

    for i in range(len(finalized_vcf.data)):
        fields = finalized_vcf.data[i]
        finalized_vcf.data[i] = fields[:9]
        finalized_vcf.data[i][8] = "GT"

    # Find the largest number of ATL alleles observed from a given locus.
    n = 0
    for fields in finalized_vcf.data:
        alt = fields[4].split(",")
        if len(alt) > n:
            n = len(alt)

    # Create fake samples in the VCF data.
    for i in range(n):
        finalized_vcf.header.append(f"TEST_SAMPLE{i + 1}")

    for i in range(len(finalized_vcf.data)):
        fields = finalized_vcf.data[i]
        v = parse_vcf_fields(fields)
        sep = "|"
        if len(v["alt"]) > 1:
            for j in range(n):
                if j + 1 > len(v["alt"]):
                    finalized_vcf.data[i].append(f"0{sep}1")
                else:
                    finalized_vcf.data[i].append(f"0{sep}{j + 1}")
        else:
            for j in range(n):
                finalized_vcf.data[i].append(f"0{sep}1")

    samples = vcf2samples(finalized_vcf, filter=False)
    snp_list = []
    
    for name in samples:
        snp_list += samples[name].hap[0].obs
        snp_list += samples[name].hap[1].obs
    
    # remove duplicates
    snp_list = list(set(snp_list))
    
    # remove non-variants
    snp_list = [x for x in snp_list if x.wt != x.var]

    # get candidates
    cand_list = [v for k, v in stardb.items() if set(v.core).issubset(snp_list) and not v.sv]

    temp = []

    temp.append(['name', 'score', 'core', 'tag', 'callable'])

    for name, star in stardb.items():

        if star.core:
            core = ",".join([x.summary() for x in star.core])
        else:
            core = "."

        if star.tag:
            tag = ",".join([x.summary() for x in star.tag])
        else:
            tag = "."

        fields = [name, str(star.score), core, tag]

        if star in cand_list:
            fields.append("yes")
        else:
            fields.append("no")

        temp.append(fields)

    result = ""

    for fields in temp:
        result += "\t".join(fields) + "\n"

    return result