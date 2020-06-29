from typing import Optional, List

from .common import VCFFile

def merge(vcfs: List[str], region: Optional[str] = None) -> str:
    """
    Merge VCF files.
    
    Returns:
        str: String representation of VCF file.

    Args:
        vcfs (str): VCF filenames.
        region (str, optional): Genomic region (chr:start:end).
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
