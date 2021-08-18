import pkgutil
from io import BytesIO

import pandas as pd
import numpy as np

class UnsupportedGeneError(Exception):
    """Raise if specified gene is not present in the gene table."""

def _check_gene(gene):
    df = load_gene_table()
    if gene not in df.Gene.unique():
        raise UnsupportedGeneError(gene)

def get_activity_score(gene, allele):
    """
    Return predicted activity score.

    The method can handle alleles with structural variation including gene
    deletion, duplication, and tandem arrangement.

    The method will return ``NaN`` for alleles with uncertain function as
    well as for alleles from a gene that does not use the activity score
    system.

    Parameters
    ----------
    gene : str
        Gene name.
    allele : str
        Star allele.

    Returns
    -------
    float
        Activity score.

    Examples
    --------

    >>> import pypgx
    >>> pypgx.get_activity_score('CYP2D6', '*1')      # Allele with normal function
    1.0
    >>> pypgx.get_activity_score('CYP2D6', '*1x2')    # Gene duplication of *1
    2.0
    >>> pypgx.get_activity_score('CYP2D6', '*1x4')    # Gene multiplication of *1
    4.0
    >>> pypgx.get_activity_score('CYP2D6', '*4')      # Allele with no function
    0.0
    >>> pypgx.get_activity_score('CYP2D6', '*4x2')    # Gene duplication of *4
    0.0
    >>> pypgx.get_activity_score('CYP2D6', '*22')     # Allele with uncertain function
    nan
    >>> pypgx.get_activity_score('CYP2D6', '*22x2')   # Gene duplication of *22
    nan
    >>> pypgx.get_activity_score('CYP2D6', '*36+*10') # Tandem arrangement
    0.25
    >>> pypgx.get_activity_score('CYP2B6', '*1')      # CYP2B6 does not have activity score
    nan
    """
    _check_gene(gene)

    df = load_activity_table()

    if gene not in df.Gene.unique():
        return np.nan

    df = df[df.Gene == gene]

    def lookup(x):
        return df[df.StarAllele == x].ActivityScore.values[0]

    def parsecnv(x):
        if 'x' in x:
            l = allele.split('x')
            base = l[0]
            times = int(l[1])
            return lookup(base) * times
        else:
            return lookup(x)

    return sum([parsecnv(x) for x in allele.split('+')])

def get_phenotype(gene, a, b):
    """
    Return predicted phenotype.

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
    _check_gene(gene)

    df = load_phenotype_table()

    if gene in df.Gene.unique():
        df = df[df.Gene == gene]
        l = [f'{a}/{b}', f'{b}/{a}']
        i = df.Diplotype.isin(l)
        return df[i].Phenotype.values[0]

    df = load_equation_table()

    if gene in df.Gene.unique():
        def one_row(r, score):
            return eval(r.Equation)
        score = get_activity_score(gene, a) + get_activity_score(gene, b)
        i = df.apply(one_row, args=(score,), axis=1)
        return df[i].Phenotype.values[0]

    raise ValueError('Something went wrong')

def load_activity_table():
    """
    Return the activity table.

    Returns
    -------
    pandas.DataFrame
        Activity table.

    Examples
    --------

    >>> import pypgx
    >>> df = pypgx.load_activity_table()
    """
    b = BytesIO(pkgutil.get_data(__name__, 'data/activity-table.csv'))
    return pd.read_csv(b)

def load_equation_table():
    """
    Return the phenotype equation table.

    Returns
    -------
    pandas.DataFrame
        Phenotype equation table.

    Examples
    --------

    >>> import pypgx
    >>> df = pypgx.load_equation_table()
    """
    b = BytesIO(pkgutil.get_data(__name__, 'data/equation-table.csv'))
    return pd.read_csv(b)

def load_gene_table():
    """
    Return the gene table.

    Returns
    -------
    pandas.DataFrame
        Activity table.

    Examples
    --------

    >>> import pypgx
    >>> df = pypgx.load_gene_table()
    """
    b = BytesIO(pkgutil.get_data(__name__, 'data/gene-table.csv'))
    return pd.read_csv(b)

def load_phenotype_table():
    """
    Return the phenotype table.

    Returns
    -------
    pandas.DataFrame
        Phenotype table.

    Examples
    --------

    >>> import pypgx
    >>> df = pypgx.load_phenotype_table()
    """
    b = BytesIO(pkgutil.get_data(__name__, 'data/phenotype-table.csv'))
    return pd.read_csv(b)
