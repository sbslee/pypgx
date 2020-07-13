from io import StringIO
from typing import List

from .bam2sdf import bam2sdf
from .sdf2gdf import sdf2gdf
from .common import sm_tag

def bam2gdf(
        gb: str,
        tg : str,
        cg: str,
        bam: List[str]
    ) -> str:
    """
    Create GDF file from BAM file(s).

    Returns:
        str: GDF file.

    Args:
        gb (str): Genome build (hg19, hg38).
        tg (str): Target gene.
        cg (str): Control gene or region.
        bam (list[str]): BAM file(s).
    """

    sdf = bam2sdf(gb, tg, cg, bam)
    sm = [sm_tag(x) for x in bam]
    result = sdf2gdf(None, sm, f=StringIO(sdf))
    return result
