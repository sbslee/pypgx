import pysam
import pandas as pd
from typing import Optional, List
from io import StringIO
from .common import bam_getter, sm_tag

COVERAGES = [1, 10, 20, 30, 40, 50, 100, 200, 300, 400, 500, 1000]

@bam_getter
def unicov(bed_file,
           bam_file: List[str],
           bam_dir: Optional[str] = None,
           bam_list: Optional[str] = None,
           **kwargs):

    input_files = kwargs["input_files"]
    names = [sm_tag(x) for x in input_files]
    s = pysam.depth("-a", "-b", bed_file, *input_files)
    f = StringIO(s)
    df = pd.read_table(f)
    df.columns = ["chr", "pos"] + names
    size = df.shape[0]

    dat = {"coverage": COVERAGES}

    for name in names:
        percs = []

        for coverage in COVERAGES:
            count = sum(df[name] >= coverage)
            perc = count / size * 100
            percs.append(perc)

        dat[name] = percs

    df = pd.DataFrame(dat)

    return df.to_string() + "\n"
