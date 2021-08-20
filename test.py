import unittest

import pypgx

df1 = pypgx.load_allele_table()
df2 = pypgx.load_diplotype_table()
df3 = pypgx.load_equation_table()
df4 = pypgx.load_gene_table()
df5 = pypgx.load_phenotype_table()

class TestPypgx(unittest.TestCase):

    def test_activity_table(self):
        self.assertEqual(pypgx.list_genes(), list(df1.Gene.unique()))

    def test_diplotype_table(self):
        self.assertEqual(len(df2.Gene.unique()), df4.PhenotypeMethod.value_counts()['Diplotype'])

    def test_equation_table(self):
        self.assertEqual(len(df3.Gene.unique()), df4.PhenotypeMethod.value_counts()['Score'])

    def test_priority_table(self):
        a = [x for x in pypgx.list_genes() if pypgx.has_phenotype(x)]
        b = list(df5.Gene.unique())
        self.assertEqual(a, b)

if __name__ == '__main__':
    unittest.main()
