import pkgutil
import pandas as pd
from io import BytesIO

def get_phenotype(gene, a, b):
    """
    Return predicted phenotype based on CPIC guidelines.

    Parameters
    ----------
    gene : str
        Gene name.
    a, b : str
        Star allele for each haplotype. The order of alleles does not matter.

    Returns
    -------
    str
        Phenotype prediction.

    Examples
    --------

    >>> import pypgx
    >>> pypgx.get_phenotype('CYP2C19', '*1', '*2')
    'Intermediate Metabolizer'
    >>> pypgx.get_phenotype('CYP2C19', '*2', '*1')
    'Intermediate Metabolizer'
    """
    df = get_phenotype_table()
    df = df[df.Gene == gene]
    l = [f'{a}/{b}', f'{b}/{a}']
    return df[df.Diplotype.isin(l)].Phenotype.values[0]

def get_phenotype_table():
    """
    Return phenotype table based on CPIC guidelines.

    Returns
    -------
    pandas.DataFrame
        Phenotype table.

    Examples
    --------

    >>> import pypgx
    >>> df = pypgx.get_phenotype_table()
    """
    b = BytesIO(pkgutil.get_data(__name__, 'data/phenotype-table.csv'))
    return pd.read_csv(b)
