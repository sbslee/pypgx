from .common import VCFFile

def minivcf(vcf_file: str,
            region: str,
            **kwargs) -> str:
    """Slice VCF file.

    Returns:
        VCF data in text.

    Args:
        vcf_file: VCF file.
        region: Target region.
    """

    vcf = VCFFile(vcf_file)
    vcf.read(region)
    return vcf.to_str()
