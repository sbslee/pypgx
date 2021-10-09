"""
The core submodule is the main suite of tools for PGx research.
"""

import pkgutil
import pathlib
from io import BytesIO

import numpy as np
import pandas as pd
from fuc import pyvcf, common

PROGRAM_PATH = pathlib.Path(__file__).parent.parent.parent.absolute()

FUNCTION_ORDER = [
    'No Function',
    'Decreased Function',
    'Possible Decreased Function',
    'Increased Function',
    'Possible Increased Function',
    'Uncertain Function',
    'Unknown Function',
    'Normal Function',
]

class AlleleNotFoundError(Exception):
    """Raise if specified allele is not present in the allele table."""

class GeneNotFoundError(Exception):
    """Raise if specified gene is not present in the gene table."""

class NotTargetGeneError(Exception):
    """Raise if specified gene is not one of the target genes."""

class PhenotypeNotFoundError(Exception):
    """Raise if specified phenotype is not present in the phenotype table."""

def build_definition_table(gene, assembly='GRCh37'):
    """
    Build the definition table of star alleles for specified gene.

    The table will only contain star alleles that are defined by SNVs and/or
    INDELs. It will not include alleles with SV (e.g. CYP2D6*5) or alleles
    with no variants (e.g. CYP2D6*2 for GRCh37 and CYP2D6*1 for GRCh38).

    Parameters
    ----------
    gene : str
        Target gene.
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.

    Returns
    -------
    fuc.pyvcf.VcfFrame
        Definition table.

    Examples
    --------

    >>> import pypgx
    >>> vf = pypgx.build_definition_table('CYP4F2')
    >>> vf.df
      CHROM       POS         ID REF ALT QUAL FILTER      INFO FORMAT *2 *3
    0    19  15990431  rs2108622   C   T    .      .  VI=V433M     GT  0  1
    1    19  16008388  rs3093105   A   C    .      .   VI=W12G     GT  1  0
    >>> vf = pypgx.build_definition_table('CYP4F2', assembly='GRCh38')
    >>> vf.df
      CHROM       POS         ID REF ALT QUAL FILTER      INFO FORMAT *2 *3
    0    19  15879621  rs2108622   C   T    .      .  VI=V433M     GT  0  1
    1    19  15897578  rs3093105   A   C    .      .   VI=W12G     GT  1  0
    """
    if not is_target_gene(gene):
        raise NotTargetGeneError(gene)

    other = 'GRCh38' if assembly == 'GRCh37' else 'GRCh37'

    df1 = load_allele_table()
    df1 = df1[df1.Gene == gene]
    variants = []
    for i, r in df1.iterrows():
        if pd.isna(r[f'{assembly}Core']):
            continue
        for variant in r[f'{assembly}Core'].split(','):
            if variant not in variants:
                variants.append(variant)
    data = {x: [] for x in pyvcf.HEADERS}

    for i, r in df1.iterrows():
        if r.SV or pd.isna(r[f'{assembly}Core']):
            continue

        data[r.StarAllele] = [
            '0' if pd.isna(r[f'{assembly}Core']) else
            '1' if x in r[f'{assembly}Core'].split(',') else
            '0' for x in variants
        ]

    df2 = load_variant_table()
    df2 = df2[df2.Gene == gene]
    for variant in variants:
        pos = int(variant.split('-')[1])
        ref = variant.split('-')[2]
        alt = variant.split('-')[3]
        s = df2[variant == df2[f'{assembly}Name']]
        data['CHROM'].append(s.Chromosome.values[0])
        data['POS'].append(pos)
        data['ID'].append(s.rsID.values[0])
        data['REF'].append(ref)
        data['ALT'].append(alt)
        data['QUAL'].append('.')
        data['FILTER'].append('.')
        data['INFO'].append(f'VI={s.Impact.values[0]}')
        data['FORMAT'].append('GT')
    meta = [
        '##fileformat=VCFv4.1',
        '##INFO=<ID=VI,Number=1,Type=String,Description="Variant impact">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    ]
    vf = pyvcf.VcfFrame.from_dict(meta, data).sort()
    return vf

def collapse_alleles(gene, alleles, assembly='GRCh37'):
    """
    Collapse redundant candidate alleles.
    """
    results = []
    for a in alleles:
        result = False
        for b in alleles:
            if a == b:
                continue
            v1 = set(list_variants(gene, alleles=a, assembly=assembly, mode='core'))
            v2 = set(list_variants(gene, alleles=b, assembly=assembly, mode='core'))
            if v1.issubset(v2):
                result = True
                break
        results.append(result)

    return [x for i, x in enumerate(alleles) if not results[i]]

def has_phenotype(gene):
    """
    Return True if specified gene has phenotype data.

    Parameters
    ----------
    gene : str
        Gene name.

    Returns
    -------
    bool
        Whether phenotype is supported.

    Examples
    --------

    >>> import pypgx
    >>> pypgx.has_phenotype('CYP2D6')
    True
    >>> pypgx.has_phenotype('CYP4F2')
    False
    """
    if not is_target_gene(gene):
        raise NotTargetGeneError(gene)

    df = load_gene_table()

    return gene in df[~df.PhenotypeMethod.isna()].Gene.to_list()

def has_score(gene):
    """
    Return True if specified gene has activity score.

    Parameters
    ----------
    gene : str
        Gene name.

    Returns
    -------
    bool
        Whether activity score is supported.

    Examples
    --------

    >>> import pypgx
    >>> pypgx.has_score('CYP2D6')
    True
    >>> pypgx.has_score('CYP2B6')
    False
    """
    if not is_target_gene(gene):
        raise NotTargetGeneError(gene)

    df = load_gene_table()

    return gene in df[df.PhenotypeMethod == 'Score'].Gene.unique()

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

def get_default_allele(gene, assembly='GRCh37'):
    """
    Get the default allele of specified gene.

    Parameters
    ----------
    gene : str
        Gene name.
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.

    Returns
    -------
    str
        Default allele.

    Examples
    --------

    >>> import pypgx
    >>> pypgx.get_default_allele('CYP2D6')
    '*2'
    >>> pypgx.get_default_allele('CYP2D6', assembly='GRCh38')
    '*1'
    """
    df = load_gene_table()
    allele = df[df.Gene == gene][f'{assembly}Default'].values[0]
    return allele

def get_function(gene, allele):
    """
    Get matched function from the allele table.

    Parameters
    ----------
    gene : str
        Gene name.
    allele : str
        Star allele.

    Returns
    -------
    str
        Function status.

    Examples
    --------

    >>> import pypgx
    >>> pypgx.get_function('CYP2D6', '*1')
    'Normal Function'
    >>> pypgx.get_function('CYP2D6', '*4')
    'No Function'
    >>> pypgx.get_function('CYP2D6', '*22')
    'Uncertain Function'
    >>> pypgx.get_function('UGT1A1', '*80+*37')
    'Decreased Function'
    >>> pypgx.get_function('CYP2D6', '*140')
    nan
    """
    if not is_target_gene(gene):
        raise NotTargetGeneError(gene)

    df = load_allele_table()
    df = df[(df.Gene == gene) & (df.StarAllele == allele)]

    if df.empty:
        raise AlleleNotFoundError(gene, allele)

    return df.Function.values[0]

def get_paralog(gene):
    """
    Get the paralog of specified gene.

    Parameters
    ----------
    gene : str
        Gene name.

    Returns
    -------
    str
        Paralog gene. Empty string if none exists.

    Examples
    --------

    >>> import pypgx
    >>> pypgx.get_paralog('CYP2D6')
    'CYP2D7'
    >>> pypgx.get_paralog('CYP2B6')
    ''
    """
    df = load_gene_table()
    df = df[df.Gene == gene]
    paralog = df.Paralog.values[0]
    if pd.isna(paralog):
        paralog = ''
    return paralog

def get_priority(gene, phenotype):
    """
    Get matched priority from the phenotype table.

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
    if not is_target_gene(gene):
        raise NotTargetGeneError(gene)

    if phenotype not in list_phenotypes():
        raise PhenotypeNotFoundError(phenotype)

    df = load_phenotype_table()
    i = (df.Gene == gene) & (df.Phenotype == phenotype)

    return df[i].Priority.values[0]

def get_ref_allele(gene, assembly='GRCh37'):
    """
    Get the reference allele for target gene.

    Parameters
    ----------
    gene : str
        Target gene.
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.

    Returns
    -------
    str
        Reference allele.

    Examples
    --------

    >>> import pypgx
    >>> pypgx.get_ref_allele('CYP2D6')
    '*1'
    >>> pypgx.get_ref_allele('NAT1')
    '*4'
    """
    df = load_gene_table()
    allele = df[df.Gene == gene]['RefAllele'].values[0]
    return allele

def get_region(gene, assembly='GRCh37'):
    """
    Get matched region from the gene table.

    Parameters
    ----------
    gene : str
        Gene name.
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.

    Returns
    -------
    str
        Requested region.
    """
    if gene not in list_genes(mode='all'):
        raise GeneNotFoundError(gene)

    df = load_gene_table()
    df = df[df.Gene == gene]

    return df[f'{assembly}Region'].values[0]

def get_score(gene, allele):
    """
    Get matched score from the allele table.

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

    See Also
    --------
    predict_score
        Predict activity score based on haplotype call.

    Examples
    --------

    >>> import pypgx
    >>> pypgx.get_score('CYP2D6', '*1')  # Allele with normal function
    1.0
    >>> pypgx.get_score('CYP2D6', '*4')  # Allele with no function
    0.0
    >>> pypgx.get_score('CYP2D6', '*22') # Allele with uncertain function
    nan
    >>> pypgx.get_score('CYP2B6', '*1')  # CYP2B6 does not have activity score
    nan
    """
    if not is_target_gene(gene):
        raise NotTargetGeneError(gene)

    if not has_score(gene):
        return np.nan

    df = load_allele_table()
    df = df[(df.Gene == gene) & (df.StarAllele == allele)]

    if df.empty:
        raise AlleleNotFoundError(gene + allele)

    return df.ActivityScore.values[0]

def get_variant_synonyms(gene, assembly='GRCh37'):
    """
    Get variant synonyms.

    Parameters
    ----------
    gene : str
        Target gene.
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.

    Returns
    -------
    dict
        Variant synonyms.

    Examples
    --------

    >>> import pypgx
    >>> pypgx.get_variant_synonyms('UGT1A1')
    {'2-234668879-CAT-CATAT': '2-234668879-C-CAT', '2-234668879-CAT-CATATAT': '2-234668879-C-CATAT'}
    >>> pypgx.get_variant_synonyms('CYP2D6')
    {}
    """
    df = load_variant_table()
    df = df[df.Gene == gene]
    synonyms = {}
    for i, r in df.iterrows():
        if pd.isna(r[f'{assembly}Synonym']):
            continue
        synonyms[r[f'{assembly}Synonym']] = r[f'{assembly}Name']
    return synonyms

def list_alleles(gene, variants=None, assembly='GRCh37'):
    """
    List all star alleles present in the allele table.

    Parameters
    ----------
    gene : str
        Target gene.
    variants : str or list, optional
        Only list alleles carrying specified variant(s) as a part of definition.
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.

    Returns
    -------
    list
        Requested alleles.

    Examples
    --------

    >>> import pypgx
    >>> pypgx.list_alleles('CYP4F2')
    ['*1', '*2', '*3']
    >>> pypgx.list_alleles('CYP2B6', variants=['19-41515263-A-G'], assembly='GRCh37')
    ['*4', '*6', '*7', '*13', '*19', '*20', '*26', '*34', '*36', '*37', '*38']
    """
    if not is_target_gene(gene):
        raise NotTargetGeneError(gene)

    df = load_allele_table()
    df = df[df.Gene == gene]

    if variants is not None:
        if isinstance(variants, str):
            variants = [variants]

        def one_row(r):
            l = []
            if not pd.isna(r[f'{assembly}Core']):
                l += r[f'{assembly}Core'].split(',')
            if not pd.isna(r[f'{assembly}Tag']):
                l += r[f'{assembly}Tag'].split(',')
            return all([x in l for x in variants])

        df = df[df.apply(one_row, axis=1)]

    return df.StarAllele.to_list()

def list_functions(gene=None):
    """
    List all functions present in the allele table.

    Parameters
    ----------
    gene : str, optional
        Return only functions belonging to this gene.

    Returns
    -------
    list
        Available functions.

    Examples
    --------

    >>> import pypgx
    >>> pypgx.list_functions()
    [nan, 'Normal Function', 'Uncertain Function', 'Increased Function', 'Decreased Function', 'No Function', 'Unknown Function', 'Possible Decreased Function', 'Possible Increased Function']
    >>> pypgx.list_functions(gene='CYP2D6')
    ['Normal Function', 'No Function', 'Decreased Function', 'Uncertain Function', 'Unknown Function', nan]
    """
    df = load_allele_table()

    if gene is not None:
        if not is_target_gene(gene):
            raise NotTargetGeneError(gene)

        df = df[df.Gene == gene]

    return list(df.Function.unique())

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

def list_phenotypes(gene=None):
    """
    List all phenotypes present in the phenotype table.

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
    df = load_phenotype_table()

    if gene is not None:
        if not is_target_gene(gene):
            raise NotTargetGeneError(gene)
        df = df[df.Gene == gene]

    return sorted(list(df.Phenotype.unique()))

def list_variants(gene, alleles=None, mode='all', assembly='GRCh37'):
    """
    List variants that are used to define star alleles.

    Some alleles, such as reference alleles, may return an empty list because
    they do not contain any variants.

    Parameters
    ----------
    gene : str
        Target gene.
    alleles : str or list, optional
        Allele name or list of alleles.
    mode : {'all', 'core', 'tag'}, default: 'all'
        Whether to return all variants, core variants only, or tag variants
        only.
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.

    Returns
    -------
    list
        Coordinate sorted list of variants.

    Examples
    --------

    >>> import pypgx
    >>> pypgx.list_variants('CYP4F2')
    ['19-15990431-C-T', '19-16008388-A-C']
    >>> pypgx.list_variants('CYP4F2', alleles=['*2'])
    ['19-16008388-A-C']
    >>> pypgx.list_variants('CYP4F2', alleles=['*2', '*3'])
    ['19-15990431-C-T', '19-16008388-A-C']
    >>> pypgx.list_variants('CYP4F2', alleles=['*2'], assembly='GRCh38')
    ['19-15897578-A-C']
    >>> pypgx.list_variants('CYP4F2', alleles=['*1'])
    []
    >>> pypgx.list_variants('CYP2B6', alleles=['*6'], mode='all')
    ['19-41495755-T-C', '19-41496461-T-C', '19-41512841-G-T', '19-41515263-A-G']
    >>> pypgx.list_variants('CYP2B6', alleles=['*6'], mode='core')
    ['19-41512841-G-T', '19-41515263-A-G']
    >>> pypgx.list_variants('CYP2B6', alleles=['*6'], mode='tag')
    ['19-41495755-T-C', '19-41496461-T-C']
    """
    if not is_target_gene(gene):
        raise NotTargetGeneError(gene)

    allele_table = load_allele_table()
    allele_table = allele_table[allele_table.Gene == gene]

    core_variants = []
    tag_variants = []

    if alleles is None:
        alleles = list_alleles(gene, assembly=assembly)

    if isinstance(alleles, str):
        alleles = [alleles]

    for allele in alleles:
        df = allele_table[allele_table.StarAllele == allele]

        if df.empty:
            raise AlleleNotFoundError(gene, allele)

        c = df[f'{assembly}Core'].values[0]
        t = df[f'{assembly}Tag'].values[0]

        if not pd.isna(c):
            core_variants += c.split(',')

        if not pd.isna(t):
            tag_variants += t.split(',')

    if mode == 'all':
        results = core_variants + tag_variants
    elif mode == 'core':
        results = core_variants
    elif mode == 'tag':
        results = tag_variants
    else:
        raise ValueError(f'Incorrect mode: {mode}')

    return common.sort_variants(set(results))

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

def load_cnv_table():
    """
    Load the CNV table.

    Returns
    -------
    pandas.DataFrame
        Requested table.

    Examples
    --------

    >>> import pypgx
    >>> df = pypgx.load_cnv_table()
    """
    b = BytesIO(pkgutil.get_data(__name__, 'data/cnv-table.csv'))
    return pd.read_csv(b)

def load_diplotype_table():
    """
    Load the diplotype table.

    Returns
    -------
    pandas.DataFrame
        Requested table.

    Examples
    --------

    >>> import pypgx
    >>> df = pypgx.load_diplotype_table()
    """
    b = BytesIO(pkgutil.get_data(__name__, 'data/diplotype-table.csv'))
    return pd.read_csv(b)

def load_equation_table():
    """
    Load the phenotype equation table.

    Returns
    -------
    pandas.DataFrame
        Requested table.

    Examples
    --------

    >>> import pypgx
    >>> df = pypgx.load_equation_table()
    """
    b = BytesIO(pkgutil.get_data(__name__, 'data/equation-table.csv'))
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

def load_phenotype_table():
    """
    Load the phenotype table.

    Returns
    -------
    pandas.DataFrame
        Requested table.

    Examples
    --------

    >>> import pypgx
    >>> df = pypgx.load_phenotype_table()
    """
    b = BytesIO(pkgutil.get_data(__name__, 'data/phenotype-table.csv'))
    return pd.read_csv(b)

def load_variant_table():
    """
    Load the variant table.

    Returns
    -------
    pandas.DataFrame
        Requested table.

    Examples
    --------

    >>> import pypgx
    >>> df = pypgx.load_phenotype_table()
    """
    b = BytesIO(pkgutil.get_data(__name__, 'data/variant-table.csv'))
    df = pd.read_csv(b)
    df.Chromosome = df.Chromosome.astype(str)
    return df

def predict_phenotype(gene, a, b):
    """
    Predict phenotype based on two haplotype calls.

    The method can handle star alleles with structural variation including
    gene deletion, duplication, and tandem arrangement.

    Parameters
    ----------
    gene : str
        Target gene.
    a, b : str
        Star allele for each haplotype. The order of alleles does not matter.

    Returns
    -------
    str
        Phenotype prediction.

    Examples
    --------

    >>> import pypgx
    >>> pypgx.predict_phenotype('CYP2D6', '*4', '*5')   # Both alleles have no function
    'Poor Metabolizer'
    >>> pypgx.predict_phenotype('CYP2D6', '*5', '*4')   # The order of alleles does not matter
    'Poor Metabolizer'
    >>> pypgx.predict_phenotype('CYP2D6', '*1', '*22')  # *22 has uncertain function
    'Indeterminate'
    >>> pypgx.predict_phenotype('CYP2D6', '*1', '*1x2') # Gene duplication
    'Ultrarapid Metabolizer'
    >>> pypgx.predict_phenotype('CYP2B6', '*1', '*4')   # *4 has increased function
    'Rapid Metabolizer'
    """
    if not is_target_gene(gene):
        raise NotTargetGeneError(gene)

    df = load_gene_table()
    phenotype_method = df[df.Gene == gene].PhenotypeMethod.values[0]

    if phenotype_method == 'Score':
        df = load_equation_table()
        df = df[df.Gene == gene]
        def one_row(r, score):
            return eval(r.Equation)
        score = predict_score(gene, a) + predict_score(gene, b)
        if np.isnan(score):
            return 'Indeterminate'
        i = df.apply(one_row, args=(score,), axis=1)
        return df[i].Phenotype.values[0]
    elif phenotype_method == 'Diplotype':
        df = load_diplotype_table()
        df = df[df.Gene == gene]
        l = [f'{a}/{b}', f'{b}/{a}']
        i = df.Diplotype.isin(l)
        return df[i].Phenotype.values[0]
    else:
        return 'Indeterminate'

def predict_score(gene, allele):
    """
    Predict activity score based on haplotype call.

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

    See Also
    --------
    get_score
        Get matched data from the allele table.

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
    if not is_target_gene(gene):
        raise NotTargetGeneError(gene)

    if not has_score(gene):
        return np.nan

    df = load_allele_table()
    df = df[df.Gene == gene]

    def parsecnv(x):
        if 'x' in x:
            l = x.split('x')
            base = l[0]
            times = int(l[1])
            return get_score(gene, base) * times
        else:
            return get_score(gene, x)

    return sum([parsecnv(x) for x in allele.split('+')])

def sort_alleles(
    alleles, by='priority', gene=None, assembly='GRCh37'
):
    """
    Sort star alleles either by priority or name.

    When ``by`` is 'priority' the method will report high priority alleles
    first. This means alleles are sorted by allele function (e.g.
    'No Function' > 'Normal Function') and then by the number of core
    variants (e.g. three SNVs > one SNV). The priority of allele function
    decreases in the following order: 'No Function', 'Decreased Function',
    'Possible Decreased Function', 'Increased Function', 'Possible Increased
    Function', 'Uncertain Function', 'Unknown Function', 'Normal Function'.

    When ``by`` is 'name' the method will report alleles with a smaller
    number first. This means, for example, '\*4' will come before '\*10'
    whereas lexicographic sorting would produce the opposite result. This is
    particularly useful when forming a diplotype (e.g. '\*4/\*10' vs.
    '\*10/\*4').

    Parameters
    ----------
    alleles : list
        List of alleles.
    by : {'priority', 'name'}, default: 'priority'
        Determines which method to use for sorting alleles:

        * 'priority': Report high priority alleles first.
        * 'name': Report alleles with a smaller number first.

    gene : str
        Target gene. Required when ``method`` is 'priority'.
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly. Used when ``method`` is 'priority'.

    Returns
    -------
    list
        Sorted list of alleles.

    Examples
    --------

    Assume we have following alleles for the *CYP2D6* gene:

    >>> alleles = ['*1', '*2', '*4', '*10']

    We can sort them by their prioirty with ``method='priority'``:

    >>> import pypgx
    >>> alleles = pypgx.sort_alleles(alleles, by='priority', gene='CYP2D6', assembly='GRCh37')
    >>> alleles
    ['*4', '*10', '*1', '*2']

    We can restore the original order by sorting again with ``method='name'``:

    >>> alleles = pypgx.sort_alleles(alleles, by='name')
    >>> alleles
    ['*1', '*2', '*4', '*10']
    """
    def func1(allele):
        if gene is None:
            raise ValueError('Target gene missing')
        function = get_function(gene, allele)
        a = FUNCTION_ORDER.index(function)
        core = list_variants(gene, alleles=allele, assembly=assembly, mode='core')
        b = len(core) * -1
        c = allele == get_ref_allele(gene, assembly=assembly) or allele == get_default_allele(gene, assembly=assembly)
        return (a, b, c)

    def func2(allele):
        cn = 1
        if '*' not in allele:
            n = 999
        else:
            _ = allele.split('+')[0].split('x')[0].replace('*', '')
            if [x for x in _ if not x.isdigit()]:
                n = 999
            else:
                n = int(''.join([x for x in _ if x.isdigit()]))
            if 'x' in allele.split('+')[0]:
                cn = int(allele.split('+')[0].split('x')[1])
        return (n, cn, len(allele), allele)

    funcs = {'priority': func1, 'name': func2}

    return sorted(alleles, key=funcs[by])
