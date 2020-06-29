from .common import VCFFile

def minivcf(fn: str, region: str) -> str:
    """
    Slice VCF file.
    
    Returns:
        str: Text version of VCF file.

    Args:
        fn (str): VCF filename.
        region (str): Genomic region.
    """
    vcf = VCFFile(fn)
    vcf.read(region)
    result = vcf.to_str()
    vcf.close()
    return result