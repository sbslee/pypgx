from .common import VCFFile

def minivcf(fn: str, region: str) -> str:
    """
    Slice VCF file.

    Returns:
        str: VCF file.

    Args:
        fn (str): VCF file.
        region (str): Target region.
    """

    vcf = VCFFile(fn)
    vcf.read(region)
    result = vcf.to_str()
    vcf.close()
    return result