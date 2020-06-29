from .common import VCFFile

def minivcf(fn: str, region: str) -> str:
    """
    Slice VCF file.
    
    Returns:
        str: String representation of VCF file.

    Args:
        fn (str): VCF filename.
        region (str): Genomic region.
    """

    vcf = VCFFile(fn)
    vcf.read(region)
    result = vcf.to_str()
    vcf.close()
    return result