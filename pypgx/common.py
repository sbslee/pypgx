import logging
import io
import pkgutil
import gzip
from typing import Optional, TextIO, List, Dict
import os
import copy

import pysam

LINE_BREAK1 = "-" * 70
LINE_BREAK2 = "*" * 70

CODES = {
    "vi0": "vi_unknown_impact",
    "vi1": "vi_low_impact",
    "vi2": "vi_moderate_impact",
    "vi3": "vi_high_impact",
    "so0": "so_unknown_variant",
    "so1": "so_synonymous_variant",
    "so2": "so_missense_variant",
    "so3": "so_frameshift_variant",
    "so4": "so_intron_variant",
    "so5": "so_splice_donor_variant",
    "so6": "so_inframe_deletion",
    "so7": "so_stop_gained",
    "so8": "so_downstream_gene_variant",
    "so9": "so_3_prime_UTR_variant",
    "so10": "so_upstream_gene_variant",
    "so11": "so_splice_acceptor_variant",
    "so12": "so_inframe_insertion",
    "so13": "so_start_lost",
    "so14": "so_5_prime_UTR_variant",
    "so15": "so_stop_lost",
    "fe0": "fe_unknown_effect",
    "fe1": "fe_no_effect",
    "rv0": "rv_unknown",
    "rv1": "rv_true",
    "rv2": "rv_false",
}

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def str2file(x):
    return io.StringIO(x)

def sort_star_names(names):
    def f(x):
        cn = 1
        if '*' not in x or x == '*DEL':
            n = 999
        else:
            n = int(''.join([y for y in x.split('+')[0].split('x')[0] if y.isdigit()]))
            if 'x' in x.split('+')[0]:
                cn = int(x.split('+')[0].split('x')[1])
        return (n, cn, len(x))

    return sorted(names, key = f)

def read_gene_table():
    result = {}
    text = pkgutil.get_data(__name__, "resources/sg/gene_table.txt").decode()
    for line in text.strip().split("\n"):
        fields = line.split("\t")
        gene = fields[1]
        if gene == "name":
            header = fields
            continue
        result[gene] = dict(zip(header, fields))
    return result

def read_pt_table():
    result = {}
    text = pkgutil.get_data(__name__, "resources/sg/pt_table.txt").decode()
    for line in text.strip().split("\n"):
        fields = line.split("\t")
        gene = fields[0]
        name = fields[2]
        rules = fields[3]
        if gene == "gene":
            header = fields
            continue
        if gene not in result:
            result[gene] = {}
        result[gene][name] = dict(zip(header, fields))
    return result

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

def sm_tag(bam: str) -> str:
    """
    Extract SM tag from BAM file.

    Returns:
        str: SM tag.

    Args:
        bam (str): BAM file.
    """

    header = pysam.view("-H", bam).strip().split("\n")

    l = []

    for line in header:
        fields = line.split("\t")
        if "@RG" == fields[0]:
            for field in fields:
                if "SM:" in field:
                    l.append(field.replace("SM:", ""))

    l = list(set(l))

    if not l:
        raise ValueError(f"SM tag not found: {bam}")

    if len(l) > 1:
        logger.warning(
            f"Multiple SM tags found (will return the first one): {bam}")
        result = l[0]
    else:
        result = l[0]

    return result




def is_chr(bam: str) -> bool:
    """
    Check whether SN tags in BAM file contain "chr" string.

    Returns:
        bool: True if found.

    Args:
        bam (str): BAM file.
    """

    header = pysam.view("-H", bam).strip().split("\n")

    l = []

    for line in header:
        fields = line.split("\t")
        if "@SQ" == fields[0]:
            for field in fields:
                if "SN:" in field:
                    l.append(field.replace("SN:", ""))

    return any(["chr" in x for x in l])

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

    def __init__(self, fn: str, f: Optional[TextIO] = None) -> None:
        """
        Initialize VCF file object.

        Args:
            fn (str): VCF file.
            f (TextIO, optional): VCF file.
        """

        if fn:
            if ".gz" in fn:
                self.f = gzip.open(fn, "rt")
            else:
                self.f = open(fn)
        else:
            self.f = f
        self.meta = []
        self.header = []
        self.data = []

    def read(self, region: Optional[str] = None) -> None:
        """
        Read VCF file.

        Args:
            region (str, optional): Target region.
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
                if chr != r[0] or pos < r[1]:
                    continue
                if pos > r[2]:
                    break
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

class StarAllele:
    def __init__(self):
        self.name = ""
        self.score = -100.0
        self.core = []
        self.tag = []
        self.sv = ""

    @property
    def ranked_as(self):
        """
        When sorted, unknown function alleles should be broken ties with 
        normal function alleles using other attributes. Increased function 
        alleles should come before normal function alleles.
        """
        if self.score < 0:
            return 1.0
        elif self.score > 1:
            return 0.99
        else:
            return self.score

    @property
    def rank(self):
        return (
            self.ranked_as,
            -1 * int(bool(self.sv)),
            -1 * len(self.core)
        )

    def __str__(self):
        return self.name

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return self.name == other.name

    def __hash__(self):
        return hash(self.name)

class SNPAllele:
    def __init__(self):
        self.pos = "" # reference genome position
        self.wt = "" # wild type (*1) allele
        self.var = "" # variant allele
        self.rs = "" # rs ID
        self.het = False # heterozygous
        self.ad = 0 # allelic depth
        self.td = 0 # total depth
        self.n = "" # SNP table number
        self.hg = "" # reference genome allele
        self.so = "" # sequence ontology
        self.fe = "" # functional effect
        self.vi = "" # variant impact
        self.rv = "" # reverting variant
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

    def summary(self):
        return (
            f"<{self.pos}:{self.wt}>{self.var}:"
            f"{self.ad}/{self.td}:{self.af:.2f}:"
            f"{get_codes_value(self.so)}:{get_codes_value(self.rv)}:"
            f"{get_codes_value(self.vi)}:{get_codes_value(self.fe)}>"
        )

def get_codes_key(value):
    for k, v in CODES.items():
        if value == v:
            return k
    return value

def get_codes_value(key):
    if key in CODES:
        return CODES[key]
    return key

def read_snp_table(
        fn: Optional[str] = None
    ) -> Dict[str, Dict[str, Dict[str, str]]]:
    """
    Create SNP allele table.

    Returns:
        dict[str, dict[str, Dict[str, str]]]: SNP allele table.

    Args:
        fn (str, optional): SNP allele table file.
    """

    if fn:
        t = fn
    else:
        t = f"{os.path.dirname(__file__)}/resources/sg/snp_table.txt"

    result = {}

    with open(t) as f:
        header = next(f).strip().split("\t")
        for line in f:
            fields = line.strip().split("\t")
            gene = fields[0]
            sg_id = fields[1]
            if gene not in result:
                result[gene] = {}
            result[gene][sg_id] = dict(zip(header, fields))

    genes = read_gene_table()

    # Some genes do not have SNPs at all (SVs only).
    for gene in [k for k, v in genes.items() if v["type"] == "target"]:
        if gene not in result:
            result[gene] = {}

    return result

def build_snpdb(
        tg: str,
        fn: Optional[str] = None
    ) -> List[SNPAllele]:
    """
    Create SNP allele database.

    Returns:
        List[SNPAllele]: SNP allele database.

    Args:
        tg (str): Target gene.
    """

    snp_table = read_snp_table(fn)

    result = []

    for k, v in snp_table[tg].items():
        snp = SNPAllele()
        snp.n = k
        snp.id = v["rs_id"]
        snp.pos = v["hg19_pos"]
        snp.hg = v["hg19_allele"]
        snp.var = v["var_allele"]
        snp.wt = v["wt_allele"]
        snp.fe = get_codes_key(v["functional_effect"])
        snp.so = get_codes_key(v["sequence_ontology"])
        snp.vi = get_codes_key(v["variant_impact"])
        snp.rv = get_codes_key(v[f"hg19_revertant"])
        snp.data = v
        result.append(snp)

    return result

def read_star_table(
        fn: Optional[str] = None
    ) -> Dict[str, Dict[str, Dict[str, str]]]:
    """
    Create star allele table.

    Returns:
        dict[str, dict[str, Dict[str, str]]]: Star allele table.

    Args:
        fn (str, optional): Star allele table file.
    """

    if fn:
        t = fn
    else:
        t = f"{os.path.dirname(__file__)}/resources/sg/star_table.txt"

    result = {}

    with open(t) as f:
        header = next(f).strip().split("\t")
        for line in f:
            fields = line.strip().split("\t")
            gene = fields[0]
            name = fields[2]
            if gene not in result:
                result[gene] = {}
            result[gene][name] = dict(zip(header, fields))

    return result

def build_stardb(
        tg: str,
        fn1: Optional[str] = None,
        fn2: Optional[str] = None,
    ) -> Dict[str, StarAllele]:

    """
    Create Star allele database.

    Returns:
        dict[str, SNPAllele]: Star allele database.

    Args:
        tg (str): Target gene.
        fn1 (str): Star allele table file.
        fn2 (str): SNP allele table file.
    """

    star_table = read_star_table(fn1)
    snpdb = build_snpdb(tg, fn2)

    result = {}

    for k, v in star_table[tg].items():
        if not v[f"hg19_has"]:
            continue

        star = StarAllele()
        star.name = k
        star.score = float(v["score"])
        star.core = [] if v["hg19_core"] in ["ref", "."] else copy.deepcopy([x for x in snpdb if f"{x.pos}:{x.wt}>{x.var}" in v[f"hg19_core"].split(",")])
        star.tag = [] if v["hg19_tag"] == "." else copy.deepcopy([x for x in snpdb if f"{x.pos}:{x.wt}>{x.var}" in v[f"hg19_tag"].split(",")])
        star.sv = "" if v["sv"] == "." else v["sv"]
        result[k] = star

    return result

def vcf2samples(vcf, filter):
    samples = {}
    for name in vcf.header[9:]:
        sample = Sample()
        sample.name = name
        i = vcf.header.index(name)
        for fields in vcf.data:
            v = parse_vcf_fields(fields)

            if filter and not any(["PS=D" in x for x in v["info"]]):
                continue


            gt = [int(x) for x in fields[i].split(":")[0].split("|")]
            al = [v["ref"]] + v["alt"]
            
            vi_list = ["no_change"] + [x for x in v["info"] if "VI=" in x][0].replace("VI=", "").split(",")
            so_list = ["no_change"] + [x for x in v["info"] if "SO=" in x][0].replace("SO=", "").split(",")
            fe_list = ["no_change"] + [x for x in v["info"] if "FE=" in x][0].replace("FE=", "").split(",")
            rv_list = ["no_change"] + [x for x in v["info"] if "RV=" in x][0].replace("RV=", "").split(",")
            
            for j in [0, 1]:
                snp = SNPAllele()
                snp.pos = v["pos"]
                snp.wt = v["ref"]
                snp.var = al[gt[j]]
                snp.rs = v["id"]
                snp.het = gt[0] != gt[1]
                snp.so = so_list[gt[j]]
                snp.vi = vi_list[gt[j]]
                snp.fe = fe_list[gt[j]]
                snp.rv = rv_list[gt[j]]
                
                if "AD" in v["format"]:
                    ad_list = [int(x) for x in fields[i].split(":")[1].split(",")]
                    snp.ad = ad_list[gt[j]]
                    snp.td = sum(ad_list)
                    
                sample.hap[j].obs.append(snp)
        samples[name] = sample
    return samples

class Haplotype:
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
        filtered = [x.af for x in self.obs if x.td > 10 and x.het and x in [y for y in self.cand[0].core]]
        return -1 if not filtered else statistics.mean(filtered)

    @property
    def af_mean_gene(self):
        filtered = [x.af for x in self.obs if x.td > 10 and x.het and self.start <= int(x.pos) <= self.end]
        return -1 if not filtered else statistics.mean(filtered)

    def remove_star(self, sX):
        """Remove the given star allele from the candidates list."""
        for i, star in enumerate(self.cand):
            if star.name == sX.name:
                del self.cand[i]
                break

    def add_dup(self, cn):
        """Duplicate the main star allele by the given CN."""
        if cn == 1:
            return
        if cn > 10:
            cn = 10
        sX = self.cand[0]
        name = sX.name + 'x' + str(cn)
        sY = StarAllele()
        
        sY.name = name; sY.score = sX.score * cn; sY.core = copy.deepcopy(sX.core); sY.sv = 'cnv{}'.format(cn)
        
        self.cand.insert(0, sY)
        self.remove_star(sX)

class Sample:
    def __init__(self):
        self.name = '' # sample ID
        self.gt = False # true if genotyped
        self.sv = ['no_sv', 'no_sv'] # one SV call per haplotype
        self.pt = '' # predicted phenotype
        self.ssr = -100.0 # sum of squared residuals
        self.dip_cand = [] # candidate stars
        self.hap = [Haplotype(), Haplotype()]

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