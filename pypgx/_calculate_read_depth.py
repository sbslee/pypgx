import pysam
from pypgx.sdk import get_sn_tags, get_sm_tags, Locus
from io import StringIO
import pandas as pd
from pypgx.sdk import Results

def calculate_read_depth(target_gene,
                         control_gene,
                         bam_path,
                         genome_build="hg19",
                         output_file=None):
    """Create a GDF (GATK DepthOfCoverage Format) file for Stargazer from
    BAM files by computing read depth.

    Parameters
    ----------
    target_gene : str
        Name of the target gene. Choices: {'abcb1', 'cacna1s',
        'cftr', 'cyp1a1', 'cyp1a2', 'cyp1b1', 'cyp2a6',
        'cyp2a13', 'cyp2b6', 'cyp2c8', 'cyp2c9', 'cyp2c19',
        'cyp2d6', 'cyp2e1', 'cyp2f1', 'cyp2j2', 'cyp2r1',
        'cyp2s1', 'cyp2w1', 'cyp3a4', 'cyp3a5', 'cyp3a7',
        'cyp3a43', 'cyp4a11', 'cyp4a22', 'cyp4b1', 'cyp4f2',
        'cyp17a1', 'cyp19a1', 'cyp26a1', 'dpyd', 'g6pd',
        'gstm1', 'gstp1', 'gstt1', 'ifnl3', 'nat1', 'nat2',
        'nudt15', 'por', 'ptgis', 'ryr1', 'slc15a2',
        'slc22a2', 'slco1b1', 'slco1b3', 'slco2b1', 'sult1a1',
        'tbxas1', 'tpmt', 'ugt1a1', 'ugt1a4', 'ugt2b7',
        'ugt2b15', 'ugt2b17', 'vkorc1', 'xpc'}.
    control_gene : str
        Name of a preselected control gene. Used for
        intrasample normalization during copy number analysis
        by Stargazer. Choices: {'egfr', 'ryr1', 'vdr'}.
        Alternatively, you can provide a custom genomic region
        with the 'chr:start-end' format (e.g. chr12:48232319-48301814).
    bam_path : str
        Read BAM files from ``bam_path``, one file path per line.
    genome_build : str, default: 'hg19'
        Build of the reference genome assembly. Choices:
        {'hg19', 'hg38'}.
    output_file : str, optional
        Path to the output file.

    Returns
    -------
    Results
        Results instance which has the following attributes: ``df``.

    """

    bam_files = []
    with open(bam_path) as f:
        for line in f:
            bam_files.append(line.strip())

    sn_tags = []
    sm_tags = []

    for bam_file in bam_files:
        sn_tags += get_sn_tags(bam_file)
        _sm_tags = get_sm_tags(bam_file)
        if not _sm_tags:
            raise ValueError(f"SM tags not found: {bam_file}")
        elif len(_sm_tags) > 1:
            raise ValueError(f"Multiple SM tags ({_sm_tags}) "
                             f"found: {bam_file}")
        else:
            sm_tags.append(list(_sm_tags)[0])

    if any(["chr" in x for x in list(set(sn_tags))]):
        chr = "chr"
    else:
        chr = ""

    loci = [Locus.from_input(target_gene, genome_build),
            Locus.from_input(control_gene, genome_build)]

    depth_data = ""

    for locus in sorted(loci, key=lambda x: x.region):
        depth_data += pysam.depth("-a", "-Q", "1", "-r",
                                  f"{chr}{locus.region}", *bam_files)

    df = pd.read_csv(StringIO(depth_data), sep="\t", header=None)

    df.columns = ["chrom", "pos"] + ["Depth_for_" + x for x in sm_tags]

    df.insert(0, "Locus", df["chrom"].astype(str) + ":" + df["pos"].astype(str))

    df.drop(columns=["chrom", "pos"], inplace=True)

    df.insert(1, "Total_Depth", df.iloc[:, 1:].sum(axis=1))
    df.insert(2, "Average_Depth_sample", df.iloc[:, 2:].mean(axis=1))

    results = Results(df=df)

    if output_file:
        df.to_csv(output_file, sep='\t', index=False)

    return results
