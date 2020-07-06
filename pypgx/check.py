from .common import logging

logger = logging.getLogger(__name__)

def check(star: str, snp: str) -> None:
    """
    Check table files for Stargazer.

    Args:
        star (str): Star allele table file.
        snp (str): SNP table file.
    """
    star_data = []
    snp_data = []
    gene_list = []

    with open(star) as f:
        star_header = next(f).strip().split("\t")
        for line in f:
            fields = line.strip().split("\t")
            gene = fields[star_header.index("gene")]
            if gene not in gene_list:
                gene_list.append(gene)
            star_data.append(fields)

    with open(snp) as f:
        snp_header = next(f).strip().split("\t")
        for line in f:
            fields = line.strip().split("\t")
            snp_data.append(fields)

    snp_dict = {}
    star_dict = {}

    for gene in gene_list:
        snp_dict[gene] = [x for x in snp_data if x[snp_header.index("gene")] == gene]
        star_dict[gene] = [x for x in star_data if x[star_header.index("gene")] == gene]

    logger.info("Star allele table:")

    for fields in star_data:
        number = int(fields[star_header.index("number")])
        if number == 1:
            count = 0
        count += 1
        if number != count:
            raise ValueError(f"Unmatched allele count: {gene.upper()}{name} (#{number}) should count {count}")

    logger.info("    Allele counts: Pass")

    def check_def(name, x, gene, genome):
        snp_list = [x[snp_header.index(f"{genome}_pos")] + ":" + x[snp_header.index("wt_allele")] + ">" + x[snp_header.index("var_allele")] for x in snp_dict[gene]]
        seen_list = []    
        if x == "." or x == "ref":
            return
        for snp in x.split(","):
            pos = int(snp.split(":")[0])
            wt = snp.split(":")[1].split(">")[0]
            var = snp.split(":")[1].split(">")[1]
            if snp not in snp_list:
                raise ValueError(f"Unrecognized allele definition ({genome}): {gene.upper()}{name} has {snp}")
            if wt == var:
                raise ValueError(f"Incorrect allele definition ({genome}): {gene.upper()}{name} has {pos}-{wt}-{var}")
            if seen_list and seen_list[-1] > pos:
                raise ValueError(f"Allele definition is not coordinate sorted: {gene.upper()}{name} (#{number}) has {x}")
            seen_list.append(pos)

    genome = "hg19"
    for fields in star_data:
        gene = fields[star_header.index("gene")]
        core = fields[star_header.index(f"{genome}_core")]
        tag = fields[star_header.index(f"{genome}_tag")]
        name = fields[star_header.index("name")]
        check_def(name, core, gene, genome)
        check_def(name, tag, gene, genome)

    logger.info("    Allele definition for hg19: Pass")

    genome = "hg38"
    for fields in star_data:
        gene = fields[star_header.index("gene")]
        core = fields[star_header.index(f"{genome}_core")]
        tag = fields[star_header.index(f"{genome}_tag")]
        name = fields[star_header.index("name")]
        check_def(name, core, gene, genome)
        check_def(name, tag, gene, genome)

    logger.info("    Allele definition for hg38: Pass")

    def check_func(x, gene, has_cadd30, type, has_lof, genome):
        snp_list = [x[snp_header.index(f"{genome}_pos")] + ":" + x[snp_header.index("wt_allele")] + ">" + x[snp_header.index("var_allele")] for x in snp_dict[gene]]
        cadd_dict = dict(zip(snp_list, [float(x[snp_header.index("cadd")]) for x in snp_dict[gene]]))
        so_dict = dict(zip(snp_list, [x[snp_header.index("sequence_ontology")] for x in snp_dict[gene]]))
        effect_dict = dict(zip(snp_list, [x[snp_header.index("functional_effect")] for x in snp_dict[gene]]))
        lof_dict = dict(zip(snp_list, [x[snp_header.index("lof")] for x in snp_dict[gene]]))
        seen_list = []    
        if x == "." or x == "ref":
            return
        cadd30_seen = False
        for snp in x.split(","):
            pos = int(snp.split(":")[0])
            wt = snp.split(":")[1].split(">")[0]
            var = snp.split(":")[1].split(">")[1]
            if type == "core":
                if cadd_dict[snp] >= 30:
                    cadd30_seen = True
                    if has_cadd30 == "no":
                        raise ValueError(f"Incorrect CADD information: {gene.upper()}{name} (#{number}) is marked as not having CADD30, but it has {snp} with CADD >= 30 ({cadd_dict[snp]})")
                if has_lof == "no" and lof_dict[snp] == "yes":
                    raise ValueError(f"Incorrect LoF information: {gene.upper()}{name} (#{number}) is marked as not having LoF, but it has {snp} with LoF (Effect={effect_dict[snp]}; SO={so_dict[snp]})")
        if type == "core" and not cadd30_seen and has_cadd30 == "yes":
            raise ValueError(f"Incorrect CADD information: {gene.upper()}{name} (#{number}) is marked as having CADD30, but no such variant was found ({x}; {','.join({str(cadd_dict[y]) for y in x.split(',')})})")

    for fields in star_data:
        gene = fields[star_header.index("gene")]
        name = fields[star_header.index("name")]
        core = fields[star_header.index("hg19_core")]
        tag = fields[star_header.index("hg19_tag")]
        number = int(fields[star_header.index("number")])
        has_cadd30 = fields[star_header.index("has_cadd30")]
        has_lof = fields[star_header.index("has_lof")]
        score = float(fields[star_header.index("score")])
        check_func(core, gene, has_cadd30, "core", has_lof, "hg19")
        check_func(tag, gene, has_cadd30, "tag", has_lof, "hg19")
        if has_cadd30 == "yes" and score < 0:
            raise ValueError(f"Conflicting data: {gene.upper()}{name} (#{number}) is marked as having CADD30, but it has unknown score")
        if has_lof == "yes" and score < 0:
            raise ValueError(f"Conflicting data: {gene.upper()}{name} (#{number}) is marked as having LoF, but it has unknown score")

    logger.info("    CADD30 information: Pass")
    logger.info("    LoF information: Pass")
    
    logger.info("SNP table:")

    for fields in snp_data:
        gene = fields[snp_header.index("gene")]
        pos = fields[snp_header.index("hg19_pos")]
        hg = fields[snp_header.index("hg19_allele")]
        var = fields[snp_header.index("var_allele")]
        wt = fields[snp_header.index("wt_allele")]
        rev = fields[snp_header.index("hg19_revertant")] == "yes"
        impact = fields[snp_header.index("variant_impact")]
        cadd = float(fields[snp_header.index("cadd")])
        causal = fields[snp_header.index("causal")]
        
        name = f"{gene.upper()}-{pos}-{hg}-{var}-{wt}"
        number = int(fields[1].replace("sg", ""))
        current_pos = int(pos)
        if number == 1:
            previous_pos = current_pos
            previous_var = var
            count = 0
        count += 1
        if current_pos < previous_pos:
            raise ValueError(f"Variants are not coordinate sorted: {name} (#{number}) is after {previous_pos}")
        if current_pos == previous_pos and var < previous_var:
            raise ValueError(f"Variants are not alphabetically sorted: {name} is after {previous_var}")
        previous_pos = current_pos
        previous_var = var
        if number != count:
            raise ValueError(f"Unmatched variant count: {name} (#{number}) should count {count}")
        if rev and (hg != var or wt == hg or wt == var):
            raise ValueError(f"Conflicting data: {name} is marked as revertant")
        if cadd >= 30 and impact != "vi_high_impact":
            raise ValueError(f"Conflicting data: {name} has CADD >= 30 and {impact}")
        
        # check associated phenotype
        phenotype = ""
        for star in star_dict[gene]:
            if "," in star[star_header.index("hg19_core")]:
                continue
            if star[star_header.index("sv")] != ".":
                continue
            if f"{pos}:{wt}>{var}" in star[star_header.index("hg19_core")].split(","):
                phenotype = star[star_header.index("phenotype")]
        if not phenotype:
            phenotype = "."
        if causal != phenotype:
            raise ValueError(f"Incorrect phenotype information: {name} ({','.join(causal)}) should have {phenotype}")
        if causal not in [".", "unknown_function", "normal_function", "IV/Normal"] and impact != "vi_high_impact":
            raise ValueError(f"Incorrect impact information: {name} ({impact}; {causal}) should have high_impact")

    logger.info("    Variant counts: Pass")
    logger.info("    Coordinate sorted: Pass")
    logger.info("    Alphabetically sorted: Pass")
    logger.info("    Reverting variants: Pass")
    logger.info("    Variant impact: Pass")
    logger.info("    Variant phenotype: Pass")