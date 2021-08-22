import os

import pandas as pd
from fuc import pyvcf

def parse_pharmvar(fn):
    """
    Parse PharmVar gene data.

    Parameters
    ----------
    fn : str
        Gene data directory.
    """
    gene = os.path.basename(fn).split('-')[0]

    rs_dict = {}

    vfs = {'GRCh37': [], 'GRCh38': []}
    alleles = {}
    for i, assembly in enumerate(['GRCh37', 'GRCh38']):
        for r, d, f in os.walk(f'{fn}/{assembly}'):
            for file in f:
                if file.endswith('.tsv'):
                    df = pd.read_table(f'{r}/{file}', comment='#')
                if file.endswith('.vcf'):
                    vf = pyvcf.VcfFrame.from_file(f'{r}/{file}')
                    vfs[assembly].append(vf)
        chrom = vfs['GRCh37'][0].contigs[0]
        for j, r in df.iterrows():
            name = r['Haplotype Name'].replace(gene, '')
            if name not in alleles:
                alleles[name] = [[], []]
            if pd.isna(r['Variant Allele']):
                continue
            variant = f"{chrom}-{r['Variant Start']}-{r['Reference Allele']}-{r['Variant Allele']}"
            rs_dict[variant] = r['rsID']
            alleles[name][i].append(variant)

    variants = {'GRCh37': {}, 'GRCh38': {}}

    for name in alleles:
        for i, assembly in enumerate(['GRCh37', 'GRCh38']):
            for variant in alleles[name][i]:
                if variant not in variants[assembly]:
                    variants[assembly][variant] = []
                    if name not in variants[assembly][variant]:
                        variants[assembly][variant].append(name)

    for name in alleles:
        alleles[name] = [','.join(alleles[name][0]), ','.join(alleles[name][1])]

    df1 = pd.DataFrame(alleles).T
    df1.columns = ['GRCh37', 'GRCh38']
    df1 = df1.replace('', 'N/A')
    df1.to_csv(f'{gene}-allele-table.csv')

    def func(r):
        if len(r.REF) == len(r.ALT) == 1:
            return ','.join(variants[assembly][r.Name])
        elif len(r.REF) == 1 and len(r.ALT) > 1:
            name = f'{r.CHROM}-{r.POS}---{r.ALT[1:]}'
            return ','.join(variants[assembly][name])
        elif len(r.REF) > 1 and len(r.ALT) == 1:
            name = f'{r.CHROM}-{r.POS}-{r.REF[1:]}--'
            return ','.join(variants[assembly][name])
        else:
            raise ValueError('Something went wrong')

    def func2(r):
        if len(r.REF) == len(r.ALT) == 1:
            return rs_dict[r.Name]
        elif len(r.REF) == 1 and len(r.ALT) > 1:
            name = f'{r.CHROM}-{r.POS}---{r.ALT[1:]}'
            return rs_dict[name]
        elif len(r.REF) > 1 and len(r.ALT) == 1:
            name = f'{r.CHROM}-{r.POS}-{r.REF[1:]}--'
            return rs_dict[name]
        else:
            raise ValueError('Something went wrong')

    for assembly in ['GRCh37', 'GRCh38']:
        df2 = pyvcf.merge(vfs[assembly]).chr_prefix().df
        df2['Name'] = df2.apply(lambda r: f'{r.CHROM}-{r.POS}-{r.REF}-{r.ALT}', axis=1)
        df2['Alleles'] = df2.apply(func, axis=1)
        df2['rsID'] = df2.apply(func2, axis=1)
        df2.to_csv(f'{gene}-{assembly}.csv')
