import os
import csv
import sys
import datetime
import pkgutil
from typing import TextIO, Optional

from .common import get_target_genes

def _add_overview_section(genotype_table, pair_table, target_genes):
    table = (
        "<table>\n"
        "<tr>\n"
        "  <th style='width: 10%;'>No.</th>\n"
        "  <th style='width: 10%;'>Gene</th>\n"
        "  <th style='width: 20%;'>Genotype</th>\n"
        "  <th style='width: 10%;'>Total AS</th>\n"
        "  <th style='width: 30%;'>Predicted Phenotype</th>\n"
        "  <th style='width: 10%;'>Drugs</th>\n"
        "  <th style='width: 10%;'>Guidelines</th>\n"
        "</tr>\n"
    )
    words = ["normal", "unknown", "-"]
    for i, gene in enumerate(target_genes, 1):
        if any([x in genotype_table[gene]["phenotype"] for x in words]):
            color = "black"
        else:
            color = "red"
        genotype = genotype_table[gene]["hap1_main"] + "/" + genotype_table[gene]["hap2_main"]
        score = genotype_table[gene]["dip_score"]
        phenotype = genotype_table[gene]["phenotype"]
        drugs = len([x["Drug"] for x in pair_table if x["Gene"] == gene.upper()])
        guidelines = len([x for x in pair_table if x["Gene"] == gene.upper() and x["Guideline"] != "-"])
        table += (
            f"<tr style='color: {color};'>\n"
            f"  <td>{i}</td>\n"
            f"  <td>{gene.upper()}</td>\n"
            f"  <td>{genotype}</td>\n"
            f"  <td>{score}</td>\n"
            f"  <td>{phenotype}</td>\n"
            f"  <td>{drugs}</td>\n"
            f"  <td>{guidelines}</td>\n"
            "</tr>\n"
        )
    return table + "</table>"

def _add_genotypes_section(genotype_table, target_genes):
    table = (
        "<table>\n"
        "<tr>\n"
        "  <th style='width: 5%;'>No.</th>\n"
        "  <th style='width: 5%;'>Gene</th>\n"
        "  <th style='width: 20%;'>Star Allele</th>\n"
        "  <th style='width: 10%;'>AS</th>\n"
        "  <th style='width: 50%;'>SNVs/Indels</th>\n"
        "  <th style='width: 10%;'>SVs</th>\n"
        "</tr>\n"
    )
    for i, gene in enumerate(target_genes, 1):
        hmc1 = "<br />".join(genotype_table[gene]["hap1_main_core"].split(","))
        hmc2 = "<br />".join(genotype_table[gene]["hap2_main_core"].split(","))
        hap1_main = genotype_table[gene]["hap1_main"]
        hap1_score = genotype_table[gene]["hap1_score"]
        hap1_sv = genotype_table[gene]["hap1_sv"]
        hap2_main = genotype_table[gene]["hap2_main"]
        hap2_score = genotype_table[gene]["hap2_score"]
        hap2_sv = genotype_table[gene]["hap2_sv"]
        table += (
            "<tr>\n"
            f"  <td rowspan='2'>{i}</td>\n"
            f"  <td rowspan='2'>{gene.upper()}</td>\n"
            f"  <td>{hap1_main}</td>\n"
            f"  <td>{hap1_score}</td>\n"
            f"  <td>{hmc1}</td>\n"
            f"  <td>{hap1_sv}</td>\n"
            "</tr>\n"
            "<tr>\n"
            f"  <td>{hap2_main}</td>\n"
            f"  <td>{hap2_score}</td>\n"
            f"  <td>{hmc2}</td>\n"
            f"  <td>{hap2_sv}</td>\n"
            "</tr>\n"
        )
    return table + "</table>"

def _add_drugs_section(genotype_table, pair_table):
    table = (
        "<table>\n"
        "<tr>\n"
        "  <th style='width: 5%;'>No.</th>\n"
        "  <th style='width: 15%;'>Drug</th>\n"
        "  <th style='width: 5%;'>Gene</th>\n"
        "  <th style='width: 5%;'>Level</th>\n"
        "  <th style='width: 10%;'>FDA</th>\n"
        "  <th style='width: 60%;'>Guideline</th>\n"
        "</tr>\n"
    )
    for i, pair in enumerate(sorted(
        pair_table, key = lambda x: (x["Drug"].lower(), x["Gene"])), 1):

        gn = pair["Gene"].lower() # gene name

        # Gene-drug pairs are bolded if a genotype is available for the gene.
        if gn in genotype_table and genotype_table[gn]["hap1_main"] != "-":
            bold = "bold"
        else:
            bold = "normal"

        # Gene-drug pairs are in red if altered phenotype is predicted.
        if bold == "normal":
            color = "black"
        elif any([x in genotype_table[gn]["phenotype"]
            for x in ["normal", "unknown", "-"]]):
            color = "black"
        else:
            color = "red"

        table += (
        f"<tr style='font-weight: {bold}; color: {color};'>\n"
        f"  <td>{i}</td>\n"
        f"  <td>{pair['Drug']}</td>\n"
        f"  <td>{pair['Gene']}</td>\n"
        f"  <td>{pair['CPIC Level']}</td>\n"
        f"  <td>{pair['PGx on FDA Label']}</td>\n"
        f"  <td>{pair['Guideline']}</td>\n"
        "</tr>\n"
        )
    return table + "</table>"

def _add_recommendations_section(action_table):
    string = ""
    for chemical in action_table:
        for gene in action_table[chemical]:
            description = (
                f"{action_table[chemical][gene]['summary']} "
                f"[PharmGKB Link: {action_table[chemical][gene]['url']}]"
            )
            table = (
                f"<h3>{chemical}-{gene}</h3>\n"
                f"<p>{description}</p>\n"
                "<table>\n"
                "<tr>\n"
                "  <th style='width: 20%;'>Phenotype</th>\n"
                "  <th style='width: 80%;'>Recommendation</th>\n"
                "</tr>\n"
            )
            for phenotype in action_table[chemical][gene]["pt"]:
                table += (
                    "<tr>\n"
                    f"  <td>{phenotype}</td>\n"
                    f"  <td>{action_table[chemical][gene]['pt'][phenotype]}</td>\n"
                    "</tr>\n"
                )
            table += "</table>\n"
            string += table
    return string

def _read_action_table():
    result = {}
    p = os.path.dirname(__file__)

    with open(f"{p}/resources/pgkb/action_table.txt") as f:
        next(f)

        for line in f:
            fields = line.strip().split("\t")
            chemical = fields[0]
            gene = fields[1]
            url = fields[2]
            summary = fields[4]
            phenotype = fields[5]
            action = fields[6]

            if chemical not in result:
                result[chemical] = {}

            if gene not in result[chemical]:
                result[chemical][gene] = {
                    "summary": summary,
                    "url": url,
                    "pt": {},
                }

            if phenotype not in result[chemical][gene]["pt"]:
                result[chemical][gene]["pt"][phenotype] = action

    return result

def _read_pair_table():
    result1 = []
    p = os.path.dirname(__file__)

    with open(f"{p}/resources/cpic/cpicPairs.csv") as f:
        result2 = next(f).strip().strip('"').replace(
            "Date last updated: ", "")
        header = next(f).strip().split(",")
        for line in f:
            fields = line.replace(", ", "/").strip().split(",")
            t = [x if x else "-" for x in fields]
            t = dict(zip(header, t))
            result1.append(t)

    return result1, result2

def _read_genotype_file(f, fn, target_genes):
    result1 = {}
    result2 = ""

    if fn:
        f = open(fn)

    header = next(f).strip().split("\t")

    for line in f:
        fields = line.strip().split("\t")
        gene = fields[0]

        if not result2:
            result2 = fields[1]

        result1[gene] = dict(zip(header, fields))

    if fn:
        f.close()

    for target_gene in target_genes:
        if target_gene not in result1:
            result1[target_gene] = dict(zip(header, ["-" for x in header]))

    return result1, result2

def gt2html(gt_file: str,
            f: Optional[TextIO] = None,
            **kwargs) -> str:
    """
    Create HTML report using Stargazer data.

    Returns:
        str: HTML report.

    Args:
        gt_file (str): Genotype file.
        f (TextIO, optional): Genotype file.
    """

    target_genes = get_target_genes()
    action_table = _read_action_table()
    pair_table, cpic_date = _read_pair_table()
    genotype_table, sample_id = _read_genotype_file(f, gt_file, target_genes)

    for k, v in genotype_table.items():
        if v["hap1_main"] == "-":
            pt = "-"
        else:
            pt = phenotyper(k, v["hap1_main"], v["hap2_main"])
        genotype_table[k]["phenotype"] = pt

    assessed = [x for x in target_genes if genotype_table[x]["status"] != "-"]
    typed = [x for x in target_genes if genotype_table[x]["status"] == "g"]

    string = (
        "<!DOCTYPE html>\n"
        "<html>\n"
        "<head>\n"
        "<title>Stargazer Report</title>\n"
        "<style>\n"
        "* {\n"
        "  font-family: Arial, Helvetica, sans-serif;\n"
        "}\n"
        "table {\n"
        "  border-collapse: collapse;\n"
        "  width: 100%;\n"
        "  font-size: 80%;\n"
        "}\n"
        "th, td {\n"
        "  border: 1px solid black;\n"
        "  padding: 4px;\n"
        "}\n"
        "</style>\n"
        "</head>\n"
        "<body>\n"
        "<h1>Stargazer Report</h1>\n"
        "<p>\n"
        f"  Sample ID: {sample_id}<br />\n"
        f"  Date: {datetime.datetime.now()}<br />\n"
        f"  Genes examined: {len(assessed)}/{len(target_genes)}<br />\n"
        f"  Genotypes called: {len(typed)}/{len(assessed)}<br />\n"
        "</p>\n"
        "<h2>Introduction</h2>\n"
        "<p>\n"
        "  Thank you for choosing Stargazer! Stargazer is a \n"
        "  bioinformatiscs tool for predicting how a person's DNA \n"
        "  affects their response to hundreds of medications. \n"
        "  Stargazer does this by accurately calling star alleles \n"
        "  (haplotypes) in pharmacogenetic (PGx) genes, which are \n"
        "  defined by single-nucleotide variants (SNVs), small \n"
        "  insertion-deletions (indels), and/or large structural \n"
        "  variants (SVs). Once identified, these star alleles can \n"
        "  be translated to an activity score (AS), which is in turn \n"
        "  used to predict the person's drug response. Stargazer can \n"
        "  utilize genomic data from various sources including \n"
        "  next-generation sequencing (NGS) and single nucleotide \n"
        "  polymorphism (SNP) array. For NGS data, Stargazer supports \n"
        "  whole genome sequencing (WGS) as well as targeted \n"
        "  sequencing such as whole exome sequencing (WES). \n"
        "  For more details please visit the Stargazer website \n"
        "  (https://stargazer.gs.washington.edu/stargazerweb/).\n"
        "</p>\n"
        "<p>\n"
        "  This report includes hundreds of gene/drug pairs \n"
        "  (e.g., CYP2D6/codeine) with accompanying levels of evidence \n"
        "  for changing drug choice and dosing decisions. These pairs \n"
        "  are described by the Clinical Pharmacogenetics Implementation \n"
        "  Consortium (CPIC), which is an international assoication \n"
        "  whose primary goal is to facilitate the use of PGx tests for \n"
        "  patient care. Most importantly, CPIC provides detailed \n"
        "  guidelines for helping clinicians understand how available \n"
        "  genetic test results should be used to optimize drug \n"
        f"  therapy. As of {cpic_date}, there are \n"
        f"  {len(pair_table)} gene/drug pairs listed in the CPIC \n"
        "  website (https://cpicpgx.org). Finally, the Food and Drug \n"
        "  Administration (FDA) provides additional guidance by \n"
        "  requiring applicable PGx test information be included in \n"
        "  the drug labeling. These FDA-approved drug labels are \n"
        "  included in this report as well.\n"
        "</p>\n"
        "<p>\n"
        "  Disclaimer: This report is still very much in development. \n"
        "  Please do not use it other than for code testing. Thank you.\n"
        "</p>\n"
        "<h2>Sections</h2>\n"
        "<ul>\n"
        "  <li>Overview</li>\n"
        "  <li>Genotypes</li>\n"
        "  <li>Drugs</li>\n"
        "  <li>Recommendations</li>\n"
        "</ul>\n"
        "<p style='page-break-before: always;'>\n"
        "<h2>Overview</h2>\n"
        "<p>\n"
        "  PGx genes whose genotype leads to altered phenotype are \n"
        "  shown in <span style='color: red;'>red</span>.\n"
        "</p>\n"
        f"{_add_overview_section(genotype_table, pair_table, target_genes)}\n"
        "<p style='page-break-before: always;'>\n"
        "<h2>Genotypes</h2>\n"
        f"{_add_genotypes_section(genotype_table, target_genes)}\n"
        "<p style='page-break-before: always;'>\n"
        "<h2>Drugs</h2>\n"
        "<p>\n"
        "  Gene/drug pairs are shown in \n"
        "  <span style='font-weight: bold;'>bold</span> if genotype is \n"
        "  available, and in \n"
        "  <span style='font-weight: bold; color: red;'>red</span> \n"
        "  if altered phenotype is predicted.\n"
        "</p>\n"
        f"{_add_drugs_section(genotype_table, pair_table)}\n"
        "<p style='page-break-before: always;'>\n"
        "<h2>Recommendations</h2>\n"
        "<p>\n"
        "  Gene/drug pairs are shown in \n"
        "  <span style='font-weight: bold;'>bold</span> if genotype is \n"
        "  available, and in \n"
        "  <span style='font-weight: bold; color: red;'>red</span> \n"
        "  if altered phenotype is predicted.\n"
        "</p>\n"
        f"{_add_recommendations_section(action_table)}\n"
        "</body>\n"
        "</html>\n"
    )

    return string
