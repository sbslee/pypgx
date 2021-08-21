import unittest

import pypgx
import pandas as pd
import numpy as np

class TestPypgx(unittest.TestCase):

    def test_activity_table(self):
        df = pypgx.load_allele_table()
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
        for gene in pypgx.list_genes():
            if not pypgx.has_definition(gene):
                continue
            df1 = pypgx.load_allele_table()
            df2 = pypgx.load_variant_table()
            df1 = df1[df1.Gene == gene]
            df2 = df2[df2.Gene == gene]
            for assembly in ['GRCh37', 'GRCh38']:
                variants = []
                for i, r in df1.iterrows():
                    if pd.isna(r[assembly]):
                        continue
                    for variant in r[assembly].split(','):
                        if variant not in variants:
                            variants.append(variant)
                s = df2[f'{assembly}Name'].unique()
                diff = set(variants) ^ set(s[~pd.isna(s)])
                if diff:
                    raise ValueError(gene, assembly, diff)

if __name__ == '__main__':
    unittest.main()
