import pkgutil
from io import BytesIO

import numpy as np
import pandas as pd
from fuc import pyvcf

FUNCTION_ORDER = [
    'No Function',
    'Decreased Function',
    'Possible Decreased Function',
    'Increased Function',
    'Possible Increased Function',
    'Normal Function',
    'Uncertain Function',
    'Unknown Function',
]

class AlleleNotFoundError(Exception):
    """Raise if specified allele is not present in the allele table."""

class GeneNotFoundError(Exception):
    """Raise if specified gene is not present in the gene table."""

class PhenotypeNotFoundError(Exception):
    """Raise if specified phenotype is not present in the phenotype table."""

def collapse_alleles(gene, alleles, assembly='GRCh37'):
    """
    Collapse redundant candidate alleles.

    Parameters
    ----------
    gene : str
        Gene name.
    alleles : list
        List of allele names.
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.

    Returns
    -------
    list
        Collapsed list.

    Examples
    --------

    >>> import pypgx
    >>> pypgx.build_definition_table('CYP4F2').df
      CHROM       POS         ID REF ALT QUAL FILTER                           INFO FORMAT *1 *2 *3
    0    19  15989040     rs1272   G   C    .      .  VI=nan;SO=3 Prime UTR Variant     GT  1  1  1
    1    19  15990431  rs2108622   C   T    .      .   VI=V433M;SO=Missense Variant     GT  0  0  1
    2    19  16008388  rs3093105   A   C    .      .    VI=W12G;SO=Missense Variant     GT  0  1  0
    >>> pypgx.collapse_alleles('CYP4F2', ['*1', '*2', '*3'])
    ['*2', '*3']
    """
    results = []

    for a in alleles:
        result = False
        for b in alleles:
            if a == b:
                continue
            if is_subset(gene, a, b, assembly=assembly):
                result = True
                break
        results.append(result)

    return [x for i, x in enumerate(alleles) if not results[i]]

def build_definition_table(gene, assembly='GRCh37'):
    """
    Build the allele definition table for specified gene.

    Parameters
    ----------
    gene : str
        Gene name.
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
      CHROM       POS         ID REF ALT QUAL FILTER                          INFO FORMAT *1 *2 *3
    0    19  15990431  rs2108622   C   T    .      .  VI=V433M;SO=Missense Variant     GT  0  0  1
    1    19  16008388  rs3093105   A   C    .      .   VI=W12G;SO=Missense Variant     GT  0  1  0
    >>> vf = pypgx.build_definition_table('CYP4F2', assembly='GRCh38')
    >>> vf.df
      CHROM       POS         ID REF ALT QUAL FILTER                          INFO FORMAT *1 *2 *3
    0    19  15879621  rs2108622   C   T    .      .  VI=V433M;SO=Missense Variant     GT  0  0  1
    1    19  15897578  rs3093105   A   C    .      .   VI=W12G;SO=Missense Variant     GT  0  1  0
    """
    if gene not in list_genes():
        raise GeneNotFoundError(gene)

    other = 'GRCh38' if assembly == 'GRCh37' else 'GRCh37'

    df1 = load_allele_table()
    df1 = df1[df1.Gene == gene]
    variants = []
    for i, r in df1.iterrows():
        if pd.isna(r[assembly]):
            continue
        for variant in r[assembly].split(','):
            if variant not in variants:
                variants.append(variant)
    data = {x: [] for x in pyvcf.HEADERS}
    for i, r in df1.iterrows():
        data[r.StarAllele] = [
            '0' if pd.isna(r[assembly]) else
            '1' if x in r[assembly].split(',') else
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
        data['INFO'].append(f'VI={s.Impact.values[0]};SO={s.SO.values[0]}')
        data['FORMAT'].append('GT')
    meta = [
        '##fileformat=VCFv4.1',
        '##INFO=<ID=VI,Number=1,Type=String,Description="Variant impact">',
        '##INFO=<ID=SO,Number=1,Type=String,Description="Sequence ontology">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    ]
    vf = pyvcf.VcfFrame.from_dict(meta, data).sort()
    return vf

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
    >>> pypgx.get_default_allele('CYP3A5')
    '*3'
    >>> pypgx.get_default_allele('CYP3A5', assembly='GRCh38')
    '*1'
    >>> pypgx.get_default_allele('CYP4F2')
    ''
    """
    df = load_gene_table()
    allele = df[df.Gene == gene][f'{assembly}Default'].values[0]
    if pd.isna(allele):
        allele = ''
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
    if gene not in list_genes():
        raise GeneNotFoundError(gene)

    df = load_allele_table()
    df = df[(df.Gene == gene) & (df.StarAllele == allele)]

    if df.empty:
        raise AlleleNotFoundError(gene, allele)

    return df.Function.values[0]

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
    if gene not in list_genes():
        raise GeneNotFoundError(gene)

    if phenotype not in list_phenotypes():
        raise PhenotypeNotFoundError(phenotype)

    df = load_phenotype_table()
    i = (df.Gene == gene) & (df.Phenotype == phenotype)

    return df[i].Priority.values[0]

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
    if gene not in list_genes():
        raise GeneNotFoundError(gene)

    if not has_score(gene):
        return np.nan

    df = load_allele_table()
    df = df[(df.Gene == gene) & (df.StarAllele == allele)]

    if df.empty:
        raise AlleleNotFoundError(gene + allele)

    return df.ActivityScore.values[0]

def has_definition(gene):
    """
    Return True if specified gene has allele definition.

    Parameters
    ----------
    gene : str
        Gene name.

    Returns
    -------
    bool
        Whether allele definition is supported.

    Examples
    --------

    >>> import pypgx
    >>> pypgx.has_definition('NUDT15')
    True
    >>> pypgx.has_definition('CYP2D6')
    False
    """
    if gene not in list_genes():
        raise GeneNotFoundError(gene)

    df = load_gene_table()

    return gene in df[~df.PharmVar.isna()].Gene.to_list()

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
    if gene not in list_genes():
        raise GeneNotFoundError(gene)

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
    if gene not in list_genes():
        raise GeneNotFoundError(gene)

    df = load_gene_table()

    return gene in df[df.PhenotypeMethod == 'Score'].Gene.unique()

def is_subset(gene, a, b, assembly='GRCh37'):
    """
    Return True if one allele is a subset of the other allele.

    The method will test whether ``a`` is a subset of ``b``; therefore,
    the order of alleles matters.

    Parameters
    ----------
    gene : str
        Gene name.
    a, b : str
        Allele name.
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.

    Returns
    -------
    bool
        Test result.

    Examples
    --------

    >>> import pypgx
    >>> pypgx.list_variants('CYP4F2', '*1')
    ['19-15989040-G-C']
    >>> pypgx.list_variants('CYP4F2', '*2')
    ['19-15989040-G-C', '19-16008388-A-C']
    >>> pypgx.is_subset('CYP4F2', '*1', '*2')
    True
    >>> pypgx.is_subset('CYP4F2', '*2', '*1')
    False
    >>> pypgx.is_subset('CYP4F2', '*1', '*1')
    True
    >>> pypgx.is_subset('CYP4F2', '*2', '*2')
    True
    """
    a = set(list_variants(gene, a, assembly=assembly))
    b = set(list_variants(gene, b, assembly=assembly))
    return a.issubset(b)

def list_alleles(gene):
    """
    List all alleles present in the allele table.

    Parameters
    ----------
    gene : str
        Gene name.

    Returns
    -------
    list
        Available alleles.

    Examples
    --------

    >>> import pypgx
    >>> pypgx.list_alleles('CYP4F2')
    ['*1', '*2', '*3']
    """
    if gene not in list_genes():
        raise GeneNotFoundError(gene)

    df = load_allele_table()
    df = df[df.Gene == gene]

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
        if gene not in list_genes():
            raise GeneNotFoundError(gene)

        df = df[df.Gene == gene]

    return list(df.Function.unique())

def list_genes():
    """
    List all genes present in the gene table.

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
        if gene not in list_genes():
            raise GeneNotFoundError(gene)
        df = df[df.Gene == gene]

    return sorted(list(df.Phenotype.unique()))

def list_variants(gene, allele, assembly='GRCh37'):
    """
    List all variants that define specified allele.

    Some alleles will return an empty list because they do not contain any
    variants (e.g. reference allele).

    Parameters
    ----------
    gene : str
        Gene name.
    allele : str
        Star allele.
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.

    Returns
    -------
    list
        Defining variants.

    Examples
    --------

    >>> import pypgx
    >>> pypgx.list_variants('CYP4F2', '*2')
    ['19-16008388-A-C']
    >>> pypgx.list_variants('CYP4F2', '*2', assembly='GRCh38')
    ['19-15897578-A-C']
    >>> pypgx.list_variants('CYP4F2', '*1')
    []
    """
    if gene not in list_genes():
        raise GeneNotFoundError(gene)

    df = load_allele_table()
    df = df[(df.Gene == gene) & (df.StarAllele == allele)]

    if df.empty:
        raise AlleleNotFoundError(gene, allele)

    variants = df[assembly].values[0]

    if pd.isna(variants):
        return []

    return variants.split(',')

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

def predict_alleles(vcf, gene, assembly='GRCh37'):
    """
    Predict candidate star alleles based on observed variants.

    The input VCF must be fully phased.

    Parameters
    ----------
    vcf : pyvcf.VcfFrame or str
        VcfFrame object or VCF file.
    gene : str
        Gene name.
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.

    Returns
    -------
    dict
        Dictionary where sample names are keys and candidate lists are values.

    Examples
    --------

    >>> from fuc import pyvcf
    >>> import pypgx
    >>> data = {
    ...     'CHROM': ['19', '19'],
    ...     'POS': [15990431, 16008388],
    ...     'ID': ['.', '.'],
    ...     'REF': ['C', 'A'],
    ...     'ALT': ['T', 'C'],
    ...     'QUAL': ['.', '.'],
    ...     'FILTER': ['.', '.'],
    ...     'INFO': ['.', '.'],
    ...     'FORMAT': ['GT', 'GT'],
    ...     'A': ['0|0', '0|1'],
    ...     'B': ['1|0', '0|1'],
    ... }
    >>> vf = pyvcf.VcfFrame.from_dict([], data)
    >>> vf.df
      CHROM       POS ID REF ALT QUAL FILTER INFO FORMAT    A    B
    0    19  15990431  .   C   T    .      .    .     GT  0|0  1|0
    1    19  16008388  .   A   C    .      .    .     GT  0|1  0|1
    >>> pypgx.predict_alleles(vf, 'CYP4F2')
    {'A': [['*1'], ['*2']], 'B': [['*3'], ['*2']]}
    """
    if gene not in list_genes():
        raise GeneNotFoundError(gene)

    table = build_definition_table(gene, assembly=assembly)

    if isinstance(vcf, str):
        vf = pyvcf.VcfFrame.from_file(vcf)
    else:
        vf = vcf

    vf = vf.filter_vcf(table)

    stars = {}

    func = lambda r: f'{r.CHROM}-{r.POS}-{r.REF}-{r.ALT}'

    for star in table.samples:
        df = table.df[table.df[star] == '1']
        if df.empty:
            stars[star] = set()
        else:
            stars[star] = set(df.apply(func, axis=1))

    samples = {}

    for sample in vf.samples:
        samples[sample] = [[], []]
        df = vf.df[sample].str.split('|', expand=True)
        df.index = vf.df.apply(func, axis=1)
        for i in [0, 1]:
            try:
                s = set(df[i][df[i] == '1'].index)
            except KeyError:
                s = set()
            for star, variants in stars.items():
                if variants.issubset(s):
                    samples[sample][i].append(star)
            samples[sample][i] = collapse_alleles(gene, samples[sample][i])
            if not samples[sample][i]:
                default = get_default_allele(gene, assembly)
                if default:
                    samples[sample][i].append(default)

    return samples

def predict_phenotype(gene, a, b):
    """
    Predict phenotype based on two haplotype calls.

    The method can handle star alleles with structural variation including
    gene deletion, duplication, and tandem arrangement.

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
    >>> pypgx.predict_phenotype('CYP2D6', '*5', '*4')   # Both alleles have no function
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
    if gene not in list_genes():
        raise GeneNotFoundError(gene)

    if has_score(gene):
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
    if gene not in list_genes():
        raise GeneNotFoundError(gene)

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

def sort_alleles(gene, alleles):
    """
    Sort candidate alleles.

    Parameters
    ----------
    gene : str
        Gene name.
    allele : str
        Star allele.

    Returns
    -------
    list
        Sorted list.
    """
    def func(x):
        x = get_function(gene, x)
        if pd.isna(x):
            i = 99
        else:
            i = FUNCTION_ORDER.index(x)
        return i
    return sorted(alleles, key=func)
