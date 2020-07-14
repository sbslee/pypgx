import gzip
import statistics
from typing import List, Dict, TextIO, Optional
from copy import deepcopy

class VCFFile:
    """
    VCF file ojbect.

    Attributes:
       f (TextIO): VCF file.
       meta (list[str]): Meta data.
       header (list[str]): Header.
       data (list[str]): Genotype data.

    Examples:

        >>> vcf = VCFFile("in.vcf") # also works with "in.vcf.gz"
        >>> vcf.read("chr10:96519437-96615962") # read CYP2C19 region only
        >>> vcf.unphase()
        >>> result = vcf.to_str()
        >>> vcf.to_file("out.vcf")
        >>> vcf.close()
    """

    def __init__(
        self,
        fn: Optional[str] = None,
        f: Optional[TextIO] = None
    ) -> None:
        """
        Initialize VCF file object.

        Args:
            fn (str, optional): VCF file.
            f (TextIO, optional): VCF file.
        """

        if fn:
            if ".gz" in fn:
                self.f = gzip.open(fn, "rt")
            else:
                self.f = open(fn)
        elif f:
            self.f = f
        else:
            self.f = None
        self.meta = []
        self.header = []
        self.data = []

    def read(
        self,
        region: Optional[str] = None,
        tidy = False,
    ) -> None:
        """
        Read VCF file.

        Args:
            region (str, optional): Target region.
            tidy (bool): Remove "chr" string.
        """

        if region:
            r = parse_region(region)
            for line in self.f:
                if "##" in line:
                    self.meta.append(line)
                    continue

                fields = line.strip().split("\t")

                if fields[0] == "#CHROM":
                    self.header = fields
                    continue

                chr = fields[0]
                pos = int(fields[1])

                if chr.replace("chr", "") != r[0].replace("chr", ""):
                    continue

                if pos < r[1]:
                    continue

                if pos > r[2]:
                    break

                if tidy:
                    fields[0] = fields[0].replace("chr", "")

                self.data.append(fields)

        else:
            for line in self.f:
                if "##" in line:
                    self.meta.append(line)
                    continue

                fields = line.strip().split("\t")

                if fields[0] == "#CHROM":
                    self.header = fields
                    continue

                if tidy:
                    fields[0] = fields[0].replace("chr", "")

                self.data.append(fields)

    def to_str(self) -> str:
        """
        Return VCF file.

        Returns:
            str: VCF file.
        """

        string = ""
        for line in self.meta:
            string += line
        string += "\t".join(self.header) + "\n"
        for fields in self.data:
            string += "\t".join(fields) + "\n"
        return string

    def to_file(self, fn: str) -> None:
        """
        Write VCF file.

        Args:
            fn (str): VCF file.
        """

        string = self.to_str()
        with open(fn, "w") as f:
            f.write(string)

    def unphase(self) -> None:
        """
        Change genotype separator from '|' to '/'.
        """

        for i in range(len(self.data)):
            self.data[i][9:] = [x.replace("|", "/") for x in self.data[i][9:]]

    def phase(self) -> None:
        """
        Change genotype separator from '/' to '|'.

        .. warning::
            This is not statistcal phasing.
        """

        for i in range(len(self.data)):
            self.data[i][9:] = [x.replace("/", "|") for x in self.data[i][9:]]

    def close(self):
        """
        Close VCF file.
        """

        self.f.close()
        self.f = None

    def copy(self, l):
        """
        Copy VCF file.
        """

        new = VCFFile()

        if "meta" in l:
            new.meta = deepcopy(self.meta)

        if "header" in l:
            new.header = deepcopy(self.header)

        if "data" in l:
            new.data = deepcopy(self.data)

        return new

    @property
    def tg(self) -> str:
        """
        Extract target gene from meta data.
        """

        result = ""

        for line in self.meta:
            if "##target_gene" in line:
                result = line.strip().replace("##target_gene=", "")
                break

        return result

    @property
    def gb(self) -> str:
        """
        Extract genome build from meta data.
        """

        result = ""

        for line in self.meta:
            if "##genome_build" in line:
                result = line.strip().replace("##genome_build=", "")
                break

        return result

class SNPAllele:
    def __init__(self):
        self.pos = '' # reference genome position
        self.wt = '' # wild type (*1) allele
        self.var = '' # variant allele
        self.rs = '' # rs ID
        self.het = False # heterozygous
        self.ad = 0 # allelic depth
        self.td = 0 # total depth
        self.n = '' # SNP table number
        self.hg = '' # reference genome allele
        self.so = '' # sequence ontology
        self.fe = '' # functional effect
        self.vi = '' # variant impact
        self.rv = '' # reverting variant
        self.gb = '' # genome build
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
        self.pt = "" # predicted phenotype
        self.ssr = -100.0 # sum of squared residuals
        self.dip_cand = [] # candidate stars
        self.hap = [BioHaplotype(), BioHaplotype()]

def vcf2biosamples(
        vcf: VCFFile,
        filter: bool = False
    ) -> List[BioSample]:
    """
    Convert a VCFFile object to a list of BioSample objects.

    Returns:
        list[BioSample]: A list of BioSample objects.

    Args:
        vcf (VCFFile): A VCFFile object.
        filter (bool): Exclude any unphased markers.
    """

    result = []

    for name in vcf.header[9:]:
        biosample = BioSample(name)
        i = vcf.header.index(name)

        for fields in vcf.data:
            v = parse_vcf_fields(fields)

            if filter and not any(["PS=D" in x for x in v["info"]]):
                continue

            gt = [int(x) for x in fields[i].split(":")[0].split("|")]
            alleles = [v["ref"]] + v["alt"]

            vi = ["NA"] + [x for x in v["info"]
                if "VI=" in x][0].replace("VI=", "").split(",")

            so = ["NA"] + [x for x in v["info"]
                if "SO=" in x][0].replace("SO=", "").split(",")

            fe = ["NA"] + [x for x in v["info"]
                if "FE=" in x][0].replace("FE=", "").split(",")

            rv = ["NA"] + [x for x in v["info"]
                if "RV=" in x][0].replace("RV=", "").split(",")

            for j in [0, 1]:
                k = gt[j]
                snpallele = SNPAllele()
                snpallele.pos = v["pos"]
                snpallele.wt = v["ref"]
                snpallele.var = alleles[k]
                snpallele.rs = v["id"]
                snpallele.het = gt[0] != gt[1]
                snpallele.so = so[k]
                snpallele.vi = vi[k]
                snpallele.fe = fe[k]
                snpallele.rv = rv[k]
                snpallele.gb = vcf.gb

                if "AD" in v["format"]:
                    ad = [int(x) for x in fields[i].split(":")[1].split(",")]
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

def parse_region(region: str) -> List[str]:
    """
    Parse region.

    Returns:
        list[str]: Parsed region [chr, start, end].

    Args:
        region (str): Region.
    """

    return (
        region.split(":")[0],
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
            chr = 23
        elif "Y" in r[0]:
            chr = 24
        else:
            chr = int(r[0].replace("chr", ""))
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
        snpallele.id = v['rs_id']
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

def read_phenotype_table(
        fn: Optional[str],
    ) -> Dict[str, Dict[str, Dict[str, str]]]:
    """
    Read phenotype table file.

    Returns:
        dict[str, dict[str, dict[str, str]]]: Phenotype table object.

    Args:
        fn (str): Phenotype table file.

    Examples:
        >>> phenotype_table = read_phenotype_table("phenotype_table.txt")
        >>> print(phenotype_table["cyp2d6"]["poor_metabolizer"]["rules"])
        ==:0,
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