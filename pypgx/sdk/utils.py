import os

import pandas as pd
from fuc import pyvcf

def parse_pharmvar(fn):
    gene = os.path.basename(fn).split('-')[0]
    alleles = {}
    vf = None
    for i, assembly in enumerate(['GRCh37', 'GRCh38']):
        for r, d, f in os.walk(f'{fn}/{assembly}'):
            for file in f:
                if file.endswith('.tsv'):
                    df = pd.read_table(f'{r}/{file}', comment='#')
                if file.endswith('.vcf') and vf is None:
                    vf = pyvcf.VcfFrame.from_file(f'{r}/{file}')
        chrom = vf.contigs[0]
        df = df[~df['Haplotype Name'].str.contains('.', regex=False)]
        for j, r in df.iterrows():
            name = r['Haplotype Name'].replace(gene, '')
            if name not in alleles:
                alleles[name] = [[], []]
            if pd.isna(r['Variant Allele']):
                continue
            variant = f"{chrom}-{r['Variant Start']}-{r['Reference Allele']}-{r['Variant Allele']}"
            alleles[name][i].append(variant)
    for name in alleles:
        alleles[name] = [','.join(alleles[name][0]), ','.join(alleles[name][1])]

    df = pd.DataFrame(alleles).T
    df.columns = ['GRCh37', 'GRCh38']
    df = df.replace('', 'N/A')
    df.to_csv(f'{gene}-allele-table.csv')
    
