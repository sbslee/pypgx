import unittest

import pypgx
import pandas as pd
import numpy as np
from fuc import pyvcf, common

class TestPypgx(unittest.TestCase):

    def test_allele_table(self):
        df = pypgx.load_allele_table()

        for assembly in ['GRCh37', 'GRCh38']:
            # Find duplicate alleles with the same definition.
            i = df[[f'{assembly}Core', 'SV']].dropna().duplicated(keep=False)
            l = sorted(df[[f'{assembly}Core', 'SV']].dropna()[i][f'{assembly}Core'].to_list())
            if l:
                raise ValueError(f'Found duplicate alleles in {assembly}: {l}')

            def one_row(r):
                if pd.isna(r[f'{assembly}Core']):
                    pass
                else:
                    ordered = ','.join(sorted(r[f'{assembly}Core'].split(','), key=lambda x: common.parse_variant(x)[1]))
                    if r[f'{assembly}Core'] != ordered:
                        raise ValueError(f"Variant are not sorted: '{r[f'{assembly}Core']}' should be changed to '{ordered}'")

            df.apply(one_row, axis=1)

        self.assertEqual(pypgx.list_genes(), list(df.Gene.unique()))

    def test_diplotype_table(self):
        df1 = pypgx.load_diplotype_table()
        df2 = pypgx.load_gene_table()
        self.assertEqual(len(df1.Gene.unique()), df2.PhenotypeMethod.value_counts()['Diplotype'])

    def test_equation_table(self):
        df1 = pypgx.load_equation_table()
        df2 = pypgx.load_gene_table()
        self.assertEqual(len(df1.Gene.unique()), df2.PhenotypeMethod.value_counts()['Score'])

    def test_priority_table(self):
        df = pypgx.load_phenotype_table()
        a = [x for x in pypgx.list_genes() if pypgx.has_phenotype(x)]
        b = list(df.Gene.unique())
        self.assertEqual(a, b)

    def test_definition_table(self):
        df1 = pypgx.load_allele_table()
        df2 = pypgx.load_variant_table()
        def one_row(r):
            for assembly in ['GRCh37', 'GRCh38']:
                other = 'GRCh38' if assembly == 'GRCh37' else 'GRCh37'
                variant = r[f'{assembly}Name']
                if pd.isna(variant):
                    return
                chrom, pos, ref, alt = common.parse_variant(variant)
                if not (chrom == r.Chromosome and
                        pos == r[f'{assembly}Position'] and
                        ref == r[f'{assembly}Allele'] and
                        (alt == r.Variant or alt == r[f'{other}Allele'])):
                        raise ValueError(f'Incorrect variant data: {variant}')
        df2.apply(one_row, axis=1)
        for gene in pypgx.list_genes():
            temp1 = df1[df1.Gene == gene]
            temp2 = df2[df2.Gene == gene]
            for assembly in ['GRCh37', 'GRCh38']:
                variants = []
                for i, r in temp1.iterrows():
                    if not pd.isna(r[f'{assembly}Core']):
                        for variant in r[f'{assembly}Core'].split(','):
                            if variant not in variants:
                                variants.append(variant)
                    if not pd.isna(r[f'{assembly}Tag']):
                        for variant in r[f'{assembly}Tag'].split(','):
                            if variant not in variants:
                                variants.append(variant)
                s = temp2[f'{assembly}Name'].unique()
                diff = set(variants) ^ set(s[~pd.isna(s)])
                if diff:
                    raise ValueError(gene, assembly, diff)

    def test_predict_alleles(self):
        a = pypgx.predict_alleles('test-data/CYP4F2-GRCh37.zip')
        b = pypgx.predict_alleles('test-data/CYP4F2-GRCh38.zip')
        self.assertEqual(['*1;', '*2;', ';', '*2:19-16008388-A-C:0.5;*1:default;'], a.data.loc['A'].to_list(), b.data.loc['A'].to_list())

if __name__ == '__main__':
    unittest.main()
