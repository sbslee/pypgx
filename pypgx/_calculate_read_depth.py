import pysam
from pypgx.sdk import get_sn_tags, get_sm_tags, Locus
from io import StringIO
import pandas as pd

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
        Name of the target gene. Choices: see ``pypgx.target_genes``.
    control_gene : str
        Name of a preselected control gene. Used for intrasample
        normalization during copy number analysis by Stargazer.
        Choices: see ``pypgx.control_genes``. Alternatively,
        you can provide a custom genomic region with the
        'chr:start-end' format (e.g. chr12:48232319-48301814)
    bam_path : str
        Read BAM files from ``bam_path``, one file path per line.
    genome_build : str, default: 'hg19'
        Build of the reference genome assembly. Choices: {'hg19', 'hg38'}.
    output_file : str, optional
        Path to the output file.

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
