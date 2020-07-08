from typing import Optional, List

from .sglib import VCFFile

def merge(
        vcfs: List[str],
        region: Optional[str] = None
    ) -> str:
    """
    Merge VCF files.

    Returns:
        str: VCF file.

    Args:
        vcfs (list[str]): VCF files.
        region (str, optional): Target region.
    """

    l = []
    for fn in vcfs:
        vcf = VCFFile(fn)
        vcf.read(region)
        l.append(vcf)
    merged = l[0]
    size = len(merged.data)

    for vcf in l[1:]:
        id = vcf.header[9]
        merged.header.append(id)
        for i in range(size):
            merged.data[i].append(vcf.data[i][9])
        vcf.close()
    
    merged.close()
    result = merged.to_str()
    return result
