from typing import List

from .bam2sdf import bam2sdf
from .sdf2gdf import sdf2gdf
from .common import sm_tag, str2file

def bam2gdf(tg : str, cg: str, bam: List[str]) -> str:
    """
    Create GDF file from BAM file(s).

    Returns:
        str: GDF file.

    Args:
        tg (str): Target gene.
        cg (str): Control gene.
        bam (list[str]): BAM file(s).
    """

    sdf = bam2sdf(tg, cg, bam)
    sm = [sm_tag(x) for x in bam]
    result = sdf2gdf(None, sm, f=str2file(sdf))
    return result