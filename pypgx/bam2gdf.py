import os
from io import StringIO
from typing import List, Optional

from .bam2sdf import bam2sdf
from .sdf2gdf import sdf2gdf
from .common import sm_tag, bam_getter

@bam_getter
def bam2gdf(
        genome_build: str,
        target_gene : str,
        control_gene: str,
        output_file: str,
        bam_file: List[str],
        bam_dir: Optional[str] = None,
        bam_list: Optional[str] = None,
        **kwargs
    ) -> None:
    """Convert BAM files to a GDF file.

    Args:
        genome_build (str):
            Genome build ('hg19' or 'hg38').
        target_gene (str):
            Name of target gene (e.g. 'cyp2d6').
        control_gene (str):
            Name or region of control gene (e.g. ‘vdr’, 
            ‘chr12:48232319-48301814’)
        output_file (str):
            Write output to this file.
        bam_file (list[str]):
            Input BAM files.
        bam_dir (str, optional):
            Use all BAM files in this directory as input.
        bam_list (str, optional):
            List of input BAM files, one file per line.
    """
    # Parse keyward arguments from the decorator.
    input_files = kwargs["input_files"]

    sdf = bam2sdf(genome_build, target_gene, control_gene, input_files)
    sm = [sm_tag(x) for x in input_files]
    result = sdf2gdf(None, sm, f=StringIO(sdf))
    with open(output_file, "w") as f:
        f.write(result)
