import copy

from .common import read_gene_table

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

class Star:
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

class SNP:
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

def _read_star_table(fn):
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

def _read_snp_table(fn):
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

    genes = read_gene_table()

    # Some genes do not have SNPs at all (SVs only).
    for gene in [k for k, v in genes.items() if v["type"] == "target"]:
        if gene not in result:
            result[gene] = {}

    return result

def liftover(star: str, snp: str, tg: str) -> str:
    """
    Convert variants in SNP table from hg19 to hg38.

    Returns:
        str: Result file.

    Args:
        star (str): Star allele table file.
        snp (str): SNP table file.
        tg (str): Target gene.
    """

    star_table = _read_star_table(star)
    snp_table = _read_snp_table(snp)

    snp_db = []
    for k, v in snp_table[tg].items():
        snp = SNP()
        snp.n = k
        snp.id = v["rs_id"]
        snp.pos = v[f"hg19_pos"]
        snp.hg = v[f"hg19_allele"]
        snp.var = v["var_allele"]
        snp.wt = v["wt_allele"]
        snp.fe = get_codes_key(v["functional_effect"])
        snp.so = get_codes_key(v["sequence_ontology"])
        snp.vi = get_codes_key(v["variant_impact"])
        snp.rv = get_codes_key(v[f"hg19_revertant"])
        snp.data = v
        snp_db.append(snp)

    # Build the star database for the target gene.
    star_db = {}
    for k, v in star_table[tg].items():
        if not v[f"hg19_has"]:
            continue
        star = Star()
        star.name = k
        star.score = float(v["score"])
        star.core = [] if v["hg19_core"] in ["ref", "."] else copy.deepcopy([x for x in snp_db if f"{x.pos}:{x.wt}>{x.var}" in v[f"hg19_core"].split(",")])
        star.tag = [] if v["hg19_tag"] == "." else copy.deepcopy([x for x in snp_db if f"{x.pos}:{x.wt}>{x.var}" in v[f"hg19_tag"].split(",")])
        star.sv = "" if v["sv"] == "." else v["sv"]
        star_db[k] = star

    result = ""

    for k, v in star_db.items():
        core = []
        for snp in v.core:
            s = "{}:{}>{}".format(snp.data["hg38_pos"], snp.data["wt_allele"], snp.data["var_allele"])
            core.append(s)
        if not core:
            core = ["."]
        tag = []
        for snp in v.tag:
            s = "{}:{}>{}".format(snp.data["hg38_pos"], snp.data["wt_allele"], snp.data["var_allele"])
            tag.append(s)
        if not tag:
            tag = ["."]
        line = ",".join(core) + "\t" + ",".join(tag)
        result += line + "\n"

    return result