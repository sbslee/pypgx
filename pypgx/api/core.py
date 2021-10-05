"""
The core submodule is the main suite of tools for PGx research.
"""

import pkgutil
from io import BytesIO

import pandas as pd

def is_target_gene(gene):
    """
    Return True if specified gene is one of the target genes.

    Parameters
    ----------
    gene : str
        Gene name.

    Returns
    -------
    bool
        True if specified gene is one of the target genes.

    Examples
    --------

    >>> import pypgx
    >>> pypgx.is_target_gene('CYP2D6')
    True
    >>> pypgx.is_target_gene('CYP2D7')
    False
    """
    return gene in list_genes(mode='target')

def list_genes(mode='target'):
    """
    List genes in the gene table.

    Parameters
    ----------
    mode : {'target', 'control', 'all'}, default: 'target'
        Specify which gene set to return.

    Returns
    -------
    list
        Gene set.

    Examples
    --------

    >>> import pypgx
    >>> pypgx.list_genes(mode='target')[:5] # First five target genes
    ['CACNA1S', 'CFTR', 'CYP1A2', 'CYP2A6', 'CYP2A13']
    >>> pypgx.list_genes(mode='control')
    ['EGFR', 'RYR1', 'VDR']
    >>> pypgx.list_genes(mode='all')[:5] # Includes pseudogenes
    ['CACNA1S', 'CFTR', 'CYP1A2', 'CYP2A6', 'CYP2A7']
    """
    df = load_gene_table()

    if mode == 'target':
        df = df[df.Target]
    elif mode == 'control':
        df = df[df.Control]
    else:
        pass

    return df.Gene.to_list()

def load_allele_table():
    """
    Load the allele table.

    Returns
    -------
    pandas.DataFrame
        Requested table.

    Examples
    --------

    >>> import pypgx
    >>> df = pypgx.load_allele_table()
    """
    b = BytesIO(pkgutil.get_data(__name__, 'data/allele-table.csv'))
    return pd.read_csv(b)

def load_gene_table():
    """
    Load the gene table.

    Returns
    -------
    pandas.DataFrame
        Requested table.

    Examples
    --------

    >>> import pypgx
    >>> df = pypgx.load_gene_table()
    """
    b = BytesIO(pkgutil.get_data(__name__, 'data/gene-table.csv'))
    return pd.read_csv(b)
