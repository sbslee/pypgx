import pkgutil
from io import BytesIO

import numpy as np
import pandas as pd

class GeneNotFoundError(Exception):
    """Raise if specified gene is not present in the gene table."""

class PhenotypeNotFoundError(Exception):
    """Raise if specified phenotype is not present in the EHR priority table."""

def predict_score(gene, allele):
    """
    Return predicted activity score.

    The method can handle star alleles with structural variation including
    gene deletion, duplication, and tandem arrangement.

    Note that the method will return ``NaN`` for alleles with uncertain
    function as well as for alleles from a gene that does not use the
    activity score system.

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
    >>> pypgx.predict_score('CYP2D6', '*1')            # Allele with normal function
    1.0
    >>> pypgx.predict_score('CYP2D6', '*1x2')          # Gene duplication of *1
    2.0
    >>> pypgx.predict_score('CYP2D6', '*1x4')          # Gene multiplication of *1
    4.0
    >>> pypgx.predict_score('CYP2D6', '*4')            # Allele with no function
    0.0
    >>> pypgx.predict_score('CYP2D6', '*4x2')          # Gene duplication of *4
    0.0
    >>> pypgx.predict_score('CYP2D6', '*22')           # Allele with uncertain function
    nan
    >>> pypgx.predict_score('CYP2D6', '*22x2')         # Gene duplication of *22
    nan
    >>> pypgx.predict_score('CYP2D6', '*36+*10')       # Tandem arrangement
    0.25
    >>> pypgx.predict_score('CYP2D6', '*1x2+*4x2+*10') # Complex event
    2.25
    >>> pypgx.predict_score('CYP2B6', '*1')            # CYP2B6 does not have activity score
    nan
    """
    if gene not in list_genes():
        raise GeneNotFoundError(gene)

    if gene not in df[df.PhenotypeMethod == 'AS'].Gene.unique():
        return np.nan

    df = load_activity_table()
    df = df[df.Gene == gene]

    def lookup(x):
        return df[df.StarAllele == x].ActivityScore.values[0]

    def parsecnv(x):
        if 'x' in x:
            l = x.split('x')
            base = l[0]
            times = int(l[1])
            return lookup(base) * times
        else:
            return lookup(x)

    return sum([parsecnv(x) for x in allele.split('+')])

def predict_phenotype(gene, a, b):
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
    >>> pypgx.predict_phenotype('CYP2D6', '*5', '*4')  # Both alleles have no function
    'Poor Metabolizer'
    >>> pypgx.predict_phenotype('CYP2D6', '*5', '*4')  # The order of alleles does not matter
    'Poor Metabolizer'
    >>> pypgx.predict_phenotype('CYP2D6', '*1', '*22') # *22 has uncertain function
    'Indeterminate'
    >>> pypgx.predict_phenotype('CYP2B6', '*1', '*4')  # *4 has increased function
    'Rapid Metabolizer'
    """
    if gene not in list_genes():
        raise GeneNotFoundError(gene)

    df = load_gene_table()
    method = df[df.Gene == gene].PhenotypeMethod.values[0]

    if method == 'AS':
        df = load_equation_table()
        df = df[df.Gene == gene]
        def one_row(r, score):
            return eval(r.Equation)
        score = predict_score(gene, a) + predict_score(gene, b)
        if np.isnan(score):
            return 'Indeterminate'
        i = df.apply(one_row, args=(score,), axis=1)
        return df[i].Phenotype.values[0]
    else:
        df = load_diplotype_table()
        df = df[df.Gene == gene]
        l = [f'{a}/{b}', f'{b}/{a}']
        i = df.Diplotype.isin(l)
        return df[i].Phenotype.values[0]

def get_priority(gene, phenotype):
    """
    Return EHR priority.

    Parameters
    ----------
    gene : str
        Gene name.
    phenotype : str
        Phenotype name.

    Returns
    -------
    str
        EHR priority.

    Examples
    --------

    >>> import pypgx
    >>> pypgx.get_priority('CYP2D6', 'Normal Metabolizer')
    'Normal/Routine/Low Risk'
    >>> pypgx.get_priority('CYP2D6', 'Ultrarapid Metabolizer')
    'Abnormal/Priority/High Risk'
    >>> pypgx.get_priority('CYP3A5', 'Normal Metabolizer')
    'Abnormal/Priority/High Risk'
    >>> pypgx.get_priority('CYP3A5', 'Poor Metabolizer')
    'Normal/Routine/Low Risk'
    """
    if gene not in list_genes():
        raise GeneNotFoundError(gene)

    if phenotype not in list_phenotypes():
        raise PhenotypeNotFoundError(phenotype)

    df = load_priority_table()
    i = (df.Gene == gene) & (df.Phenotype == phenotype)

    return df[i].Priority.values[0]

def list_genes():
    """
    Return all genes present in the gene table.

    Returns
    -------
    list
        Available genes.

    Examples
    --------

    >>> import pypgx
    >>> pypgx.list_genes()
    ['CYP2B6', 'CYP2C9', 'CYP2C19', 'CYP2D6', 'CYP3A5']
    """
    df = load_gene_table()
    return df.Gene.to_list()

def list_phenotypes(gene=None):
    """
    Return all phenotypes present in the EHR priority table.

    Parameters
    ----------
    gene : str, optional
        Return only phenotypes belonging to this gene.

    Returns
    -------
    list
        Available phenotypes.

    Examples
    --------

    >>> import pypgx
    >>> pypgx.list_phenotypes()
    ['Intermediate Metabolizer', 'Normal Metabolizer', 'Poor Metabolizer', 'Rapid Metabolizer', 'Ultrarapid Metabolizer', 'Likely Intermediate Metabolizer', 'Likely Poor Metabolizer', 'Possible Intermediate Metabolizer']
    >>> pypgx.list_phenotypes(gene='CYP2D6')
    ['Ultrarapid Metabolizer', 'Normal Metabolizer', 'Intermediate Metabolizer', 'Poor Metabolizer']
    """
    df = load_priority_table()

    if gene is not None:
        if gene not in list_genes():
            raise GeneNotFoundError(gene)
        df = df[df.Gene == gene]

    return list(df.Phenotype.unique())

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

def load_diplotype_table():
    """
    Return the diplotype table.

    Returns
    -------
    pandas.DataFrame
        diplotype table.

    Examples
    --------

    >>> import pypgx
    >>> df = pypgx.load_diplotype_table()
    """
    b = BytesIO(pkgutil.get_data(__name__, 'data/diplotype-table.csv'))
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

def load_priority_table():
    """
    Return the EHR priority table.

    Returns
    -------
    pandas.DataFrame
        EHR priority table.

    Examples
    --------

    >>> import pypgx
    >>> df = pypgx.load_priority_table()
    """
    b = BytesIO(pkgutil.get_data(__name__, 'data/priority-table.csv'))
    return pd.read_csv(b)
