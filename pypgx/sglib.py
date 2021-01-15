# This script is shared between Stargazer and PyPGx.
# Author: Seung-been "Steven" Lee
# Email: sbstevenlee@gmail.com
# Last updated: 2020-08-06 14:04

import gzip
import pandas as pd
import statistics
from typing import List, Dict, TextIO, Optional
from copy import deepcopy

class SNPAllele:
    """SNP allele object.

    Attributes:
       pos (str): Genomic coordinate.
       wt (str): Wild-type allele.
       var (str): Variant allele compared to the wild type.
       rs (str): rs number.
       het (bool): True if heterozygous.
       ad (int): Allelic depth.
       td (int): Total read depth.
       n (str): SNP table number.
       hg (str): Reference allele in the assembly.
       so (str): Sequence ontology.
       fe (str): Function effect.
       vi (str): Variant impact.
       rv (str): Reverting variant.
       gb (str): Genome build.
    """
    def __init__(self):
        self.pos = ''
        self.wt = ''
        self.var = ''
        self.rs = ''
        self.het = False
        self.ad = 0
        self.td = 0
        self.n = ''
        self.hg = ''
        self.so = ''
        self.fe = ''
        self.vi = ''
        self.rv = ''
        self.gb = ''
        self.data = {}

    @property
    def key(self):
        return (self.pos, self.wt, self.var)

    @property
    def af(self):
        return 0 if self.td == 0 else self.ad / self.td

    def __eq__(self, other):
        return self.key == other.key

    def __hash__(self):
        return hash(self.key)

    def summary(self) -> str:
        """
        Return SNP summary.

        Returns:
            str: Summary of this SNPAllele object.
        """

        return (
            f"[{self.rs}:{self.gb}:{self.pos}:{self.wt}>{self.var}:"
            f"{self.ad}/{self.td}:{self.af:.2f}:"
            f"{self.so}:{self.rv}:{self.vi}:{self.fe}]"
        )

class StarAllele:
    def __init__(self):
        self.name = ''
        self.score = -100.0
        self.core = []
        self.tag = []
        self.sv = ''

    @property
    def ranked_as(self):
        """
        Unknown function alleles should be broken ties with normal function
        alleles using attributes other than activity score. Increased
        function alleles should come before normal function alleles.
        """

        if self.score < 0:
            return 1.0
        elif self.score > 1:
            return 0.99
        else:
            return self.score

    @property
    def rank(self):
        return (self.ranked_as, -1 * int(bool(self.sv)), -1 * len(self.core))

    def __str__(self):
        return self.name

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return self.name == other.name

    def __hash__(self):
        return hash(self.name)

class BioHaplotype:
    def __init__(self):
        self.cand = []
        self.obs = []
        self.start = -1 # target gene start
        self.end = -1 # target gene end

    @property
    def sv(self):
        sv = 'no_sv'
        sv_list = []
        for star_allele in self.cand:
            if star_allele.sv and star_allele.sv not in sv_list:
                sv_list.append(star_allele.sv)
        if len(sv_list) > 1:
            raise ValueError('haplotype contains multiple SV calls')
        if len(sv_list) == 1:
            sv = sv_list[0]
        return sv

    @property
    def af(self):
        return [0] if not self.obs else [x.af for x in self.obs]

    @property
    def af_mean_main(self):
        filtered = [x.af for x in self.obs
            if x.td > 10 and x.het and x in [y for y in self.cand[0].core]]

        return -1 if not filtered else statistics.mean(filtered)

    @property
    def af_mean_gene(self):
        filtered = [x.af for x in self.obs
            if x.td > 10 and x.het and self.start <= int(x.pos) <= self.end]

        return -1 if not filtered else statistics.mean(filtered)

    def remove_star(self, sX):
        """Remove the given star allele from the candidates list."""
        for i, star in enumerate(self.cand):
            if star.name == sX.name:
                del self.cand[i]
                break

    def add_dup(self, cn):
        """
        Duplicate the main star allele by the given CN.
        """
        if cn == 1:
            return

        if cn > 10:
            cn = 10

        sX = self.cand[0]
        name = sX.name + "x" + str(cn)
        sY = StarAllele()
        sY.name = name
        sY.score = sX.score * cn
        sY.core = deepcopy(sX.core)
        sY.sv = "cnv{}".format(cn)

        self.cand.insert(0, sY)
        self.remove_star(sX)

class BioSample:
    def __init__(self, name: str) -> None:
        """
        Initialize a BioSample object.

        Args:
            name (str): Sample ID.
        """

        self.name = name
        self.gt = False # true if genotyped
        self.sv = ["no_sv", "no_sv"] # one SV call per haplotype
        self.ssr = -100.0 # sum of squared residuals
        self.dip_cand = [] # candidate stars
        self.hap = [BioHaplotype(), BioHaplotype()]
        self.ssr_df = pd.DataFrame()
        self.af_df = pd.DataFrame()

def vcf2biosamples(vcf,
                   filter: bool = False) -> List[BioSample]:
    """Convert a VCFFile to a list of BioSample.

    Returns:
        A list of BioSample.

    Args:
        vcf: A VCFFile object.
        filter: If true, exclude any unphased markers.
    """

    result = []

    for name in vcf.header[9:]:
        biosample = BioSample(name)
        i = vcf.header.index(name)

        for v in vcf.data:
            if filter and "D" not in v.info["PS"]:
                continue

            gt = [int(x) for x in v.fields[i].split(":")[0].split("|")]
            alleles = [v.ref] + v.alt
            vi = ["NA"] + v.info["VI"].split(",")
            so = ["NA"] + v.info["SO"].split(",")
            fe = ["NA"] + v.info["FE"].split(",")
            rv = ["NA"] + v.info["RV"].split(",")

            for j in [0, 1]:
                k = gt[j]
                snpallele = SNPAllele()
                snpallele.pos = str(v.pos)
                snpallele.wt = v.ref
                snpallele.var = alleles[k]
                snpallele.rs = v.id
                snpallele.het = gt[0] != gt[1]
                snpallele.so = so[k]
                snpallele.vi = vi[k]
                snpallele.fe = fe[k]
                snpallele.rv = rv[k]
                snpallele.gb = vcf.search_meta("genome_build")

                if "AD" in v.format:
                    ad = [int(x) for x in v.fields[i].split(":")[1].split(",")]
                    snpallele.ad = ad[k]
                    snpallele.td = sum(ad)

                biosample.hap[j].obs.append(snpallele)

        result.append(biosample)

    return result

def sort_star_names(names: List[str]) -> List[str]:
    """
    Sort star names.

    Returns:
        list[str]: Sorted star names.

    Args:
        names (list[str]): Star names.

    Examples:

        >>> print(sort_star_names(["*33", "*4"]))
        ['*4', '*33']
        >>> print(sort_star_names(["*3x2", "*3"]))
        ['*3', '*3x2']
        >>> print(sort_star_names(["*3x2", "*4"]))
        ['*3x2', '*4']
        >>> print(sort_star_names(["*3+*5", "*4"]))
        ['*3+*5', '*4']
    """

    def f(name):
        cn = 1

        if "*" not in name or name == "*DEL":
            n = 999

        else:
            _ = name.split("+")[0].split("x")[0]
            n = int("".join([x for x in _ if x.isdigit()]))

            if "x" in name.split("+")[0]:
                cn = int(name.split("+")[0].split("x")[1])

        return (n, cn, len(name))

    return sorted(names, key = f)

def parse_region(
        region: str,
        omit: bool = False
    ) -> List[str]:
    """
    Parse region.

    Returns:
        list[str]: Parsed region [chr, start, end].

    Args:
        region (str): Region to be parsed.
        omit (bool): Remove the 'chr' string.
    """
    if omit:
        chr = region.split(":")[0].replace("chr", "")
    else:
        chr = region.split(":")[0]

    return (
        chr,
        int(region.split(":")[1].split("-")[0]),
        int(region.split(":")[1].split("-")[1]),
    )

def sort_regions(regions: List[str]) -> List[str]:
    """
    Sort regions.

    Returns:
        list[str]: Sorted regions.

    Args:
        regions (list[str]): Regions.
    """

    def f(x):
        r = parse_region(x)
        if "X" in r[0]:
            contig = 23
        elif "Y" in r[0]:
            contig = 24
        else:
            _ = r[0].replace("chr", "")

            if _.isdigit():
                contig = int(_)
            else:
                contig = 25

        return (chr, r[1], r[2])

    return sorted(regions, key = f)

def parse_vcf_fields(fields):
    return {
        "chrom": fields[0].replace("chr", ""),
        "pos": fields[1],
        "id": fields[2],
        "ref": fields[3],
        "alt": fields[4].split(","),
        "qual": fields[5],
        "filter": fields[6].split(";"),
        "info": fields[7].split(";"),
        "format": fields[8].split(":")
    }

def read_sv_table(fn):
    result = {}

    with open(fn) as f:
        header = next(f).strip().split()
        for line in f:
            fields = line.strip().split()
            gene = fields[0]
            name = fields[2]

            if gene not in result:
                result[gene] = {}

            result[gene][name] = dict(zip(header, fields))

    return result

def read_gene_table(
        fn: str
    ) -> Dict[str, Dict[str, str]]:
    """
    Read gene table file.

    Returns:
        dict[str, dict[str, str]]: Gene table object.

    Args:
        fn (str): Gene table file.

    Examples:

        >>> gene_table = read_gene_table("gene_table.txt")
        >>> print(gene_table["cyp2d6"]["hg19_start"])
        42522500
    """

    result = {}

    with open(fn) as f:
        header = next(f).strip().split("\t")

        for line in f:
            fields = line.strip().split("\t")
            name = fields[1]
            result[name] = dict(zip(header, fields))

    return result

def read_snp_table(
        fn: str,
        gene_table: Dict[str, Dict[str, str]]
    ) -> Dict[str, Dict[str, Dict[str, str]]]:
    """
    Read SNP table file.

    Returns:
        dict[str, dict[str, dict[str, str]]]: SNP table object.

    Args:
        fn (str): SNP table file.
        gene_table (dict[str, dict[str, str]]): Gene table object.

    Examples:

        >>> gene_table = read_gene_table("gene_table.txt")
        >>> snp_table = read_snp_table("snp_table.txt", gene_table)
        >>> print(snp_table["cyp2d6"]["sg2"]["wt_allele"])
        C
    """

    result = {}

    with open(fn) as f:
        header = next(f).strip().split("\t")

        for line in f:
            fields = line.strip().split("\t")
            gene = fields[0]
            sg_id = fields[1]

            if gene not in result:
                result[gene] = {}

            result[gene][sg_id] = dict(zip(header, fields))

    target_genes = [
        x for x in gene_table if gene_table[x]["type"] == "target"]

    for target_gene in target_genes:
        if target_gene not in result:
            result[target_gene] = {}

    return result

def read_star_table(
        fn: str
    ) -> Dict[str, Dict[str, Dict[str, str]]]:
    """
    Read star table file.

    Returns:
        Dict[str, Dict[str, Dict[str, str]]]: Star table object.

    Args:
        fn (str): Star table file.

    Examples:

        >>> star_table = read_star_table("star_table.txt")
        >>> print(star_table["cyp2d6"]["*2"]["hg19_core"])
        42522613:C>G,42523943:G>A
    """

    result = {}

    with open(fn) as f:
        header = next(f).strip().split("\t")
        for line in f:
            fields = line.strip().split("\t")
            gene = fields[0]
            name = fields[2]
            if gene not in result:
                result[gene] = {}
            result[gene][name] = dict(zip(header, fields))

    return result

def build_snpdb(
        tg: str,
        gb: str,
        snp_table: Dict[str, Dict[str, Dict[str, str]]],
    ) -> List[SNPAllele]:
    """
    Build SNPAllele database for target gene.

    Returns:
        list[SNPAllele]: SNPAllele objects.

    Args:
        tg (str): Target gene.
        gb (str): Genome build.
        snp_table (dict[str, dict[str, dict[str, str]]]): SNP table object.
    """

    result = []

    for k, v in snp_table[tg].items():
        snpallele = SNPAllele()
        snpallele.n = k
        snpallele.rs = v['rs_id']
        snpallele.pos = v[f'{gb}_pos']
        snpallele.hg = v[f'{gb}_allele']
        snpallele.var = v['var_allele']
        snpallele.wt = v['wt_allele']
        snpallele.fe = v['functional_effect']
        snpallele.so = v['sequence_ontology']
        snpallele.vi = v['variant_impact']
        snpallele.rv = v[f'{gb}_revertant']
        snpallele.gb = gb
        snpallele.data = v
        result.append(snpallele)

    return result

def build_stardb(
        tg: str,
        gb: str,
        star_table: Dict[str, Dict[str, Dict[str, str]]],
        snpdb: List[SNPAllele],
    ) -> Dict[str, StarAllele]:
    """
    Build StarAllele database for target gene.

    Returns:
        dict[str, StarAllele]: StarAllele objects.

    Args:
        tg (str): Target gene.
        gb (str): Genome build.
        star_table (dict[str, dict[str, dict[str, str]]]): Star table object.
        snpdb (list[SNPAllele]): SNPAllele objects.

    Examples:

        >>> gene_table = read_gene_table("gene_table.txt")
        >>> snp_table = read_snp_table("snp_table.txt", gene_table)
        >>> snpdb = build_snpdb("cyp2d6", "hg19", snp_table)
        >>> star_table = read_star_table("star_table.txt")
        >>> stardb = build_stardb("cyp2d6", "hg19", star_table, snpdb)
        >>> print(stardb["*2"].score)
        1.0
    """

    result = {}

    for k, v in star_table[tg].items():
        if v[f"{gb}_has"] == "no":
            continue

        starallele = StarAllele()
        starallele.name = k
        starallele.score = float(v["score"])

        if v[f"{gb}_core"] in ["ref", "."]:
            starallele.core = []
        else:
            starallele.core = deepcopy([
                x for x in snpdb
                if f"{x.pos}:{x.wt}>{x.var}"  in v[f"{gb}_core"].split(",")
            ])

        if v[f"{gb}_tag"] == ".":
            starallele.tag = []
        else:
            starallele.tag = deepcopy([
                x for x in snpdb
                if f"{x.pos}:{x.wt}>{x.var}" in v[f"{gb}_tag"].split(",")
            ])

        if v["sv"] == ".":
            starallele.sv = ""
        else:
            starallele.sv = v["sv"]

        result[k] = starallele

    return result
