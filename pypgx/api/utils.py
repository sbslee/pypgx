from . import data

def get_phenotype(gene, a, b):
    df = data.get_phenotype_table()
    df = df[df.Gene == gene]
    print(df[df.Diplotype.isin([f'{a}/{b}', f'{b}/{a}'])].Phenotype.values[0])
