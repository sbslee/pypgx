from .common import VCFFile
from .liftover import SNP, _read_star_table, _read_snp_table, get_codes_key, get_codes_value, Star

import os
import copy

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
        sY = Star()
        
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
                snp = SNP()
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

def peek(tg, vcf) -> str:
    """
    Find all possible star alleles from VCF file.

    Returns:
        str: Result file.

    Args:
        tg (str): Target gene.
        vcf (str): VCF file.
    """

    # Remove sample data from the VCF file.
    finalized_vcf = VCFFile(vcf)
    finalized_vcf.read()
    finalized_vcf.header = finalized_vcf.header[:9]

    for i in range(len(finalized_vcf.data)):
        fields = finalized_vcf.data[i]
        finalized_vcf.data[i] = fields[:9]
        finalized_vcf.data[i][8] = "GT"

    # Find the largest number of ATL alleles observed from a given locus.
    n = 0
    for fields in finalized_vcf.data:
        alt = fields[4].split(",")
        if len(alt) > n:
            n = len(alt)

    # Create fake samples in the VCF data.
    for i in range(n):
        finalized_vcf.header.append(f"TEST_SAMPLE{i + 1}")

    for i in range(len(finalized_vcf.data)):
        fields = finalized_vcf.data[i]
        v = parse_vcf_fields(fields)
        sep = "|"
        if len(v["alt"]) > 1:
            for j in range(n):
                if j + 1 > len(v["alt"]):
                    finalized_vcf.data[i].append(f"0{sep}1")
                else:
                    finalized_vcf.data[i].append(f"0{sep}{j + 1}")
        else:
            for j in range(n):
                finalized_vcf.data[i].append(f"0{sep}1")

    samples = vcf2samples(finalized_vcf, filter=False)

    finalized_vcf.close()

    snp_list = []
    
    for name in samples:
        snp_list += samples[name].hap[0].obs
        snp_list += samples[name].hap[1].obs
    
    # remove duplicates
    snp_list = list(set(snp_list))
    
    # remove non-variants
    snp_list = [x for x in snp_list if x.wt != x.var]
    



    star_table = _read_star_table(f"{os.path.dirname(__file__)}/resources/sg/star_table.txt")
    snp_table = _read_snp_table(f"{os.path.dirname(__file__)}/resources/sg/snp_table.txt")

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


    # get candidates
    cand_list = [v for k, v in star_db.items() if set(v.core).issubset(snp_list) and not v.sv]

    temp = []

    temp.append(['name', 'score', 'core', 'tag', 'callable'])

    for name, star in star_db.items():

        if star.core:
            core = ",".join([x.summary() for x in star.core])
        else:
            core = "."

        if star.tag:
            tag = ",".join([x.summary() for x in star.tag])
        else:
            tag = "."

        fields = [name, str(star.score), core, tag]

        if star in cand_list:
            fields.append("yes")
        else:
            fields.append("no")

        temp.append(fields)

    result = ""

    for fields in temp:
        result += "\t".join(fields) + "\n"

    return result