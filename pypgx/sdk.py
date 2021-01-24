import pypgx
import pysam
import pandas as pd

class Locus:
    """Store various information for the given gene.

    Parameters
    ----------
    name : str
        Name of the gene.
    chrom : str
        Chromosome where the gene is located.
    gene_start : int
        Start position of the gene.
    gene_end : int
        End position of the gene.
    upstream : int
        Upstream.
    downstream : int
        Downstream.
    masked_starts : list
        Start positions of the masked regions, if present.
    masked_ends : list
        End positions of the masked regions, if present.
    exon_starts : list
        Start positions of the exons.
    exon_ends : list
        End positions of the exons.

    Attributes
    ----------
    name : str
        Name of the gene.
    chrom : str
        Chromosome where the gene is located.
    gene_start : int
        Start position of the gene.
    gene_end : int
        End position of the gene.
    upstream : int
        Upstream.
    downstream : int
        Downstream.
    region_start : int
        Start position of the bigger region.
    region_end : int
        End position of the bigger region.
    masked_starts : list
        Start positions of the masked regions, if present.
    masked_ends : list
        End positions of the masked regions, if present.
    region : str
        Genomic region in the `chr:start-end` format.
    region_size : int
        Size of the bigger region.
    gene_size : int
        Size of the gene.
    exon_starts : list
        Start positions of the exons.
    exon_ends : list
        End positions of the exons.

    """
    def __init__(self, name, genome_build, chrom, gene_start, gene_end,
                 upstream, downstream, masked_starts, masked_ends,
                 exon_starts, exon_ends):
        self.name = name
        self.genome_build = genome_build
        self.chrom = chrom
        self.gene_start = gene_start
        self.gene_end = gene_end
        self.gene_size = gene_end - gene_start
        self.upstream = upstream
        self.downstream = downstream
        self.region_start = gene_start - upstream
        self.region_end = gene_end + downstream
        self.region = f"{chrom}:{gene_start-upstream}-{gene_end+downstream}"
        self.region_size = (gene_end + downstream) - (gene_start - upstream)
        self.masked_starts = masked_starts
        self.masked_ends = masked_ends
        self.exon_starts = exon_starts
        self.exon_ends = exon_ends

    @classmethod
    def from_input(cls, text, genome_build):
        """Initiate with a user input."""
        if text in pypgx.gene_df["name"].to_list():
            return cls.from_gene(text, genome_build)
        else:
            return cls.from_region(text, genome_build)

    @classmethod
    def from_gene(cls, gene, genome_build):
        """Initiate with a gene name (e.g. 'cyp2d6')."""
        row = pypgx.gene_df[pypgx.gene_df["name"] == gene]
        return cls.from_row(row, genome_build)

    @classmethod
    def from_row(cls, row, genome_build):
        """Initiate with a ``gene_df`` row."""
        name = row["name"].values[0]
        chrom = row["chrom"].values[0].replace("chr", "")
        gene_start = row[f"{genome_build}_start"].values[0]
        gene_end = row[f"{genome_build}_end"].values[0]
        upstream = row["upstream"].values[0]
        downstream = row["downstream"].values[0]
        _masked_starts = row[f"{genome_build}_masked_starts"].values[0]
        _masked_ends = row[f"{genome_build}_masked_ends"].values[0]
        _exon_starts = row[f"{genome_build}_exon_starts"].values[0]
        _exon_ends = row[f"{genome_build}_exon_ends"].values[0]
        def f(x):
            return [int(y) for y in x.strip(',').split(',')]
        masked_starts = [] if pd.isna(_masked_starts) else f(_masked_starts)
        masked_ends = [] if pd.isna(_masked_ends) else f(_masked_ends)
        exon_starts = [] if pd.isna(_exon_starts) else f(_exon_starts)
        exon_ends = [] if pd.isna(_exon_ends) else f(_exon_ends)
        return cls(name, genome_build, chrom, gene_start, gene_end,
                   upstream, downstream, masked_starts, masked_ends,
                   exon_starts, exon_ends)

    @classmethod
    def from_region(cls, region, genome_build):
        """Initiate with a custom region (e.g. 'chr12:48232319-48301814')."""
        name = region
        chrom = region.split(':')[0].replace("chr", '')
        gene_start = int(region.split(':')[1].split('-')[0])
        gene_end = int(region.split(':')[1].split('-')[1])
        upstream = 0
        downstream = 0
        masked_starts = []
        masked_ends = []
        exon_starts = []
        exon_ends = []
        return cls(name, genome_build, chrom, gene_start, gene_end,
                   upstream, downstream, masked_starts, masked_ends,
                   exon_starts, exon_ends)

class Results:
    """Namespace for storing various results from pypgx."""
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

def get_sm_tags(bam_file):
    """Extract the SM tags (sample names) from a BAM file.

    Parameters
    ----------
    bam_file : str
        Path to the BAM file.

    Returns
    -------
    set
        SM tags.

    """
    header = pysam.view("-H", bam_file).strip().split("\n")
    _sm_tags = []
    for line in header:
        fields = line.split("\t")
        if "@RG" == fields[0]:
            for field in fields:
                if "SM:" in field:
                    _sm_tags.append(field.replace("SM:", ""))
    return set(_sm_tags)

def get_sn_tags(bam_file):
    """Extract the SN tags (sequence or contig names) from a BAM file.

    Parameters
    ----------
    bam_file : str
        Path to the BAM file.

    Returns
    -------
    set
        SN tags.

    """
    header = pysam.view("-H", bam_file).strip().split("\n")
    sn_tags = []
    for line in header:
        fields = line.split("\t")
        if "@SQ" == fields[0]:
            for field in fields:
                if "SN:" in field:
                    sn_tags.append(field.replace("SN:", ""))
    return set(sn_tags)

def get_activity_score(target_gene, haplotype_call):
    """Convert haplotype call to activity score.

    Parameters
    ----------
    target_gene : str
        Name of the target gene (e.g. cyp2d6).
    haplotype_call : str
        Haplotype call (e.g. \*2).

    Returns
    -------
    float
        Activity_score.

    """
    df = pypgx.star_df[pypgx.star_df["gene"] == target_gene]
    activity_score = 0
    for star_allele in haplotype_call.split('+'):
        if 'x' in star_allele:
            copies = int(star_allele.split('x')[1])
            name = star_allele.split('x')[0]
            activity_score += df[df["name"] == name][
                "score"].values[0] * copies
        else:
            activity_score += df[df["name"] == star_allele][
                "score"].values[0]
    return activity_score
