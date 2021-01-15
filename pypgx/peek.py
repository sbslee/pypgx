import os

from .common import get_stardb, VCFFile
from .sglib import vcf2biosamples


def peek(vcf_file: str,
         **kwargs) -> str:
    """Find all possible star alleles from VCF file.

    Returns:
        Summary of star allele status.

    Args:
        vcf_file: Stargazer VCF file (finalized.vcf).
    """

    # Remove sample data from the VCF file.
    finalized_vcf = VCFFile(vcf_file)
    finalized_vcf.read()
    finalized_vcf.header = finalized_vcf.header[:9]

    stardb = get_stardb(finalized_vcf.search_meta("target_gene"),
                        finalized_vcf.search_meta("genome_build"))

    for i in range(len(finalized_vcf.data)):
        record = finalized_vcf.data[i]
        record.fields = record.fields[:9]
        record.fields[8] = "GT"

    # Find the largest number of ATL alleles observed from a given locus.
    n = 0
    for record in finalized_vcf.data:
        if len(record.alt) > n:
            n = len(record.alt)

    # Create fake samples in the VCF data.
    for i in range(n):
        finalized_vcf.header.append(f"TEST_SAMPLE{i + 1}")

    for i in range(len(finalized_vcf.data)):
        record = finalized_vcf.data[i]
        sep = "|"
        if len(record.alt) > 1:
            for j in range(n):
                if j + 1 > len(record.alt):
                    record.fields.append(f"0{sep}1")
                else:
                    record.fields.append(f"0{sep}{j + 1}")
        else:
            for j in range(n):
                record.fields.append(f"0{sep}1")

    biosamples = vcf2biosamples(finalized_vcf, filter=False)
    snp_list = []

    for biosample in biosamples:
        snp_list += biosample.hap[0].obs
        snp_list += biosample.hap[1].obs

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
