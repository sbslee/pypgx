from .sglib import VCFFile

def minivcf(vcf_file: str, region: str, **kwargs) -> str:
    """
    Slice VCF file.

    Returns:
        str: VCF file.

    Args:
        vcf_file (str): VCF file.
        region (str): Target region.
    """

    vcf = VCFFile(vcf_file)
    vcf.read(region)
    result = vcf.to_str()
    vcf.close()
    return result