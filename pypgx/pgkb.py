import json
import time
import logging
import sys

from .common import LINE_BREAK1, LINE_BREAK2

import requests
import pandas as pd
from bs4 import BeautifulSoup

def _update(result, gene, chemical, table, url, summary, type):
    """
    Update the result dictionary.
    """

    if chemical not in result:
        result[chemical] = {}

    if gene not in result[chemical]:
        result[chemical][gene] = {
            "summary": summary,
            "url": url,
            "type": type,
            "phenotypes": {},
        }

    # Get the indicies of column names with the word "recommendation".
    indicies = []
    for column in list(table):
        if "recommendation" in column.lower():
            i = list(table).index(column)
            indicies.append(i)
    if not indicies:
        logging.warning(LINE_BREAK2)
        logging.warning(f"No recommendations found: {table}")
        logging.warning(LINE_BREAK2)
        return

    # Update the recommendation for each phenotype.
    phenotypes = table.iloc[:, 0].tolist()
    for i in range(len(phenotypes)):
        recommandation = []
        for j in indicies:
            title = list(table)[j]
            content = str(table.iloc[:,j].tolist()[i])
            recommandation.append(f"{title}: {content}")
        result[chemical][gene]["phenotypes"][phenotypes[i]] = " ".join(recommandation)

def pgkb(args):

    # Get the PharmGKB IDs of CPIC guidelines.
    pgkb_id = []
    base_url = "https://api.pharmgkb.org/v1/data/guideline"
    response = requests.get(f"{base_url}?source=cpic&view=min")
    json = response.json()
    for guideline in json["data"]:
        pgkb_id.append(guideline["id"])

    # Get the CPIC guidelines.
    result = {}
    for id in pgkb_id:
        logging.info(LINE_BREAK1)
        request_url = f"{base_url}/{id}?view=base"
        logging.info(f"Request URL: {request_url}")
        response = requests.get(request_url)
        logging.info(f"Response: {response}")
        json = response.json()
        frontend_url = json["data"]["@id"]
        logging.info(f"Frontend URL: {frontend_url}")

        # Get the guideline summary.
        markup = json["data"]["summaryMarkdown"]["html"]
        soup = BeautifulSoup(markup, "lxml")
        summary = soup.text.strip().replace("\n", " ")

        # Get the associated genes.
        genes = []
        for gene in json["data"]["relatedGenes"]:
            genes.append(gene["symbol"])

        # Get the associated chemicals.
        chemicals = []
        for chemical in json["data"]["relatedChemicals"]:
            chemicals.append(chemical["name"])

        # Get the associated tables.
        markup = json["data"]["textMarkdown"]["html"]
        soup = BeautifulSoup(markup, "lxml")
        while soup.find_all("sup"):
            sup_tag = soup.sup
            sup_tag.decompose()
        tables = pd.read_html(str(soup))

        logging.info(f"Genes: {genes}")
        logging.info(f"Chemicals: {chemicals}")
        
        ng = len(genes)
        nc = len(chemicals)
        nt = len(tables)

        # Case 1: This is the simpletest case where there are one gene,
        # one chemical and one table. Nothing much left to parse.
        if ng == nc == nt == 1:
            logging.info("Type: Case 1")
            _update(result, genes[0], chemicals[0], tables[0], frontend_url, summary, "case1")

        # Case 2: This is a simple case where there are one chemical,
        # one table and multiple genes. Nothing much left to parse.
        elif ng > 1 and nc == nt == 1:
            logging.info("Type: Case 2")
            for gene in genes:
                _update(result, gene, chemicals[0], tables[0], frontend_url, summary, "case2")

        # Case 3: This is a simple case where there are one gene
        # one table and multiple chemicals. Nothing much left to parse.
        elif nc > 1 and ng == nt == 1:
            logging.info("Type: Case 3")
            for chemical in chemicals:
                _update(result, genes[0], chemical, tables[0], frontend_url, summary, "case3")

        # Case 4: This is a simple case where there are one chemical and
        # identical number of genes and tables. Just need to find the correct
        # gene-table pairs
        elif nc == 1 and ng == nt:
            logging.info("Type: Case 4")
            temp = {}
            for table in tables:
                phenotypes = table.iloc[:, 0].tolist()
                seen = []
                for phenotype in phenotypes:
                    for gene in genes:
                        if gene in phenotype and gene not in seen:
                            seen.append(gene)
                if len(seen) == 1:
                    temp[seen[0]] = table
            if len(temp) != nt:
                logging.warning(LINE_BREAK2)
                logging.warning("Case 4: cannot find gene-table pairs")
                logging.warning(LINE_BREAK2)
                continue
            for i, gene in enumerate(genes):
                _update(result, gene, chemicals[0], temp[gene], frontend_url, summary, "case4")

        # Case 5: This is a complicated case where there are one chemical,
        # two genes and three tables. Need to find the correct
        # gene-table pairs.
        elif ng == 2 and nc == 1 and nt == 3:
            logging.info("Type: Case 5")
            temp = {}
            for table in tables:
                seen = []
                for gene in genes:
                    if gene in "".join(list(table)) and gene not in seen:
                        seen.append(gene)
                    if table.iloc[:, 0].str.contains(gene, regex=False).any() and gene not in seen:
                        seen.append(gene)
                if len(seen) == 1:
                    temp[seen[0]] = table
            if len(temp) != ng:
                logging.warning(LINE_BREAK2)
                logging.warning("Case 5: cannot find gene-table pairs")
                logging.warning(LINE_BREAK2)
                continue
            for i, gene in enumerate(genes):
                _update(result, gene, chemicals[0], temp[gene], frontend_url, summary, "case5")

        # Case 6: This guideline (voriconazole-CYP2C19) has separate
        # tables for adult and pediatric patients.
        elif id == "PA166161537":
            _update(result, genes[0], chemicals[0], tables[0], frontend_url, summary, "case6")

        # Case 7: This guideline (atomoxetine-CYP2D6) has separate
        # tables for adult and pediatric patients.
        elif id == "PA166181885":
            _update(result, genes[0], chemicals[0], tables[1], frontend_url, summary, "case7")

        # Case 8: This guideline has two genes (CACNA1S and RYR1), seven
        # chemicals and one table.
        elif id == "PA166180457":
            for gene in genes:
                for chemical in chemicals:
                    _update(result, gene, chemical, tables[0], frontend_url, summary, "case8")

        # Case 9: This guideline (warfarin-CYP2C9, CYP4F2, VKORC1) requires
        # manual review.
        elif id == "PA166104949":
            logging.warning(LINE_BREAK2)
            logging.warning(f"This guideline requires manual review")
            logging.warning(LINE_BREAK2)
            continue

        elif len(tables) > 1 and len(genes) == len(chemicals) == 1:
            logging.warning(LINE_BREAK2)
            logging.warning(f"More than one tables found: {tables}")
            logging.warning(LINE_BREAK2)
            continue

        else:
            logging.warning(LINE_BREAK2)
            logging.warning(f"Found {len(genes)} genes, {len(chemicals)} chemicals, and {len(tables)} tables")
            logging.warning(LINE_BREAK2)
            continue

        # Wait some time before the next request.
        time.sleep(2)
        logging.info(LINE_BREAK1)

    string = "chemical\tgene\turl\ttype\tsummary\tphenotype\trecommendation\n"
    for chemical in sorted(result):
        for gene in result[chemical]:
            url = result[chemical][gene]["url"]
            summary = result[chemical][gene]["summary"]
            type = result[chemical][gene]["type"]
            for phenotype in result[chemical][gene]["phenotypes"]:
                recommendation = result[chemical][gene]["phenotypes"][phenotype]
                string += f"{chemical}\t{gene}\t{url}\t{type}\t{summary}\t{phenotype}\t{recommendation}\n"

    if args.o:
        with open(args.o, "w") as f:
            f.write(string)
    else:
        sys.stdout.write(string)