from typing import Optional, List

from vcfgo.VCFFile import VCFFile

def mergevcf(vcf_file: List[str],
             region: Optional[str] = None,
             **kwargs) -> str:
    """Merge VCF files.

    Returns:
        Merged VCF file.

    Args:
        vcf_file: VCF files to be merged.
        region: Target region.
    """

    l = []

    for fn in vcf_file:
        vcf = VCFFile(fn)

        if region:
            vcf.read(region)
        else:
            vcf.read()

        l.append(vcf)

    merged = l[0]

    size = len(merged.data)

    for vcf in l[1:]:
        id = vcf.header[9]
        merged.header.append(id)

        for i in range(size):
            record = vcf.data[i]
            merged.data[i].fields.append(record.fields[9])

    return merged.to_str()
