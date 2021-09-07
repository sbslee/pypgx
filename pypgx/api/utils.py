import pkgutil
from io import BytesIO
import tempfile
import subprocess
import os
import pickle

from .. import sdk

import numpy as np
import pandas as pd
from fuc import pybam, pyvcf, pycov, common
from sklearn import model_selection
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import SVC
from sklearn.impute import SimpleImputer

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

###################
# Private methods #
###################

##################
# Public methods #
##################

def build_definition_table(gene, assembly='GRCh37'):
    """
    Build the definition table of star alleles for specified gene.

    The table will only contain star alleles that are defined by SNVs and/or
    INDELs. It will not include alleles with SV (e.g. CYP2D6*5) or alleles
    with no variants (e.g. CYP2D6*2 for GRCh37 and CYP2D6*1 for GRCh38).

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
      CHROM       POS         ID REF ALT QUAL FILTER                          INFO FORMAT *2 *3
    0    19  15990431  rs2108622   C   T    .      .  VI=V433M;SO=Missense Variant     GT  0  1
    1    19  16008388  rs3093105   A   C    .      .   VI=W12G;SO=Missense Variant     GT  1  0
    >>> vf = pypgx.build_definition_table('CYP4F2', assembly='GRCh38')
    >>> vf.df
      CHROM       POS         ID REF ALT QUAL FILTER                          INFO FORMAT *2 *3
    0    19  15879621  rs2108622   C   T    .      .  VI=V433M;SO=Missense Variant     GT  0  1
    1    19  15897578  rs3093105   A   C    .      .   VI=W12G;SO=Missense Variant     GT  1  0
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
        if r.SV or pd.isna(r[assembly]):
            continue

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

def compute_control_statistics(
    bam=None, fn=None, gene=None, region=None, assembly='GRCh37'
):
    """
    Compute copy number from read depth for target gene.

    Parameters
    ----------
    target : pypgx.Archive
        Archive file with the semantic type CovFrame[ReadDepth].
    control : pypgx.Archive
        Archive file with the semandtic type ControlStatistics.
    """
    bam_files = []

    if bam is None and fn is None:
        raise ValueError(
            "Either the 'bam' or 'fn' parameter must be provided.")
    elif bam is not None and fn is not None:
        raise ValueError(
            "The 'bam' and 'fn' parameters cannot be used together.")
    elif bam is not None and fn is None:
        if isinstance(bam, str):
            bam_files.append(bam)
        else:
            bam_files += bam
    else:
        bam_files += common.convert_file2list(fn)

    df = load_gene_table()

    if gene is not None:
        region = df[df.Gene == gene][f'{assembly}Region'].values[0]

    if all([pybam.has_chr(x) for x in bam_files]):
        region = 'chr' + region

    cf = pycov.CovFrame.from_bam(
        bam=bam_files, region=region, zero=False
    )

    metadata = {
        'Control': gene,
        'Assembly': assembly,
        'SemanticType': 'TSV[Statistics]',
    }
    data = cf.df.iloc[:, 2:].describe()
    result = sdk.Archive(metadata, data)
    return result

def compute_copy_number(target, control, samples):
    """
    Compute copy number from read depth for target gene.

    The method will convert read depth to copy number by performing
    intra-sample normalization with control statistics.

    If the input data was generated with targeted sequencing, as opposed to
    WGS, the method will also apply inter-sample normalization using the
    population statistics. However, for best results it's recommended to
    manually specify a list of known samples without SV.

    Parameters
    ----------
    target : pypgx.Archive
        Archive file with the semantic type CovFrame[ReadDepth].
    control : pypgx.Archive
        Archive file with the semandtic type TSV[Statistcs].
    samples : list, optional
        List of known samples without SV.

    Returns
    -------
    pypgx.Archive
        Archive file with the semandtic type CovFrame[CopyNumber].
    """
    depth = sdk.Archive.from_file(target)
    statistics = sdk.Archive.from_file(control)

    # Apply intra-sample normalization.
    df = depth.data.copy_df()
    medians = statistics.data.T['50%']
    df.iloc[:, 2:] = df.iloc[:, 2:] / medians * 2

    # Apply inter-sample normalization.
    if depth.metadata['Platform'] == 'Targeted':
        if samples is None:
            medians = df.iloc[:, 2:].median(axis=1).replace(0, np.nan)
        else:
            medians = df[samples].median(axis=1).replace(0, np.nan)
        df.iloc[:, 2:] = df.iloc[:, 2:].div(medians, axis=0) * 2

    cf = pycov.CovFrame(df)
    metadata = depth.copy_metadata()
    metadata['SemanticType'] = 'CovFrame[CopyNumber]'
    metadata['Control'] = statistics.metadata['Control']

    return sdk.Archive(metadata, cf)

def compute_target_depth(
    gene, bam=None, fn=None, assembly='GRCh37'
):
    """
    Compute read depth for target gene with BAM data.

    Parameters
    ----------
    target : pypgx.Archive
        Archive file with the semantic type CovFrame[ReadDepth].
    control : pypgx.Archive
        Archive file with the semandtic type ControlStatistics.
    """
    bam_files = []

    if bam is None and fn is None:
        raise ValueError(
            "Either the 'bam' or 'fn' parameter must be provided.")
    elif bam is not None and fn is not None:
        raise ValueError(
            "The 'bam' and 'fn' parameters cannot be used together.")
    elif bam is not None and fn is None:
        if isinstance(bam, str):
            bam_files.append(bam)
        else:
            bam_files += bam
    else:
        bam_files += common.convert_file2list(fn)

    if all([pybam.has_chr(x) for x in bam_files]):
        prefix = 'chr'
    else:
        prefix = ''

    region = get_region(gene, assembly=assembly)

    data = pycov.CovFrame.from_bam(
        bam=bam_files, region=f'{prefix}{region}', zero=True
    )

    metadata = {
        'Gene': gene,
        'Assembly': assembly,
        'SemanticType': 'CovFrame[ReadDepth]',
    }
    result = sdk.Archive(metadata, data)
    return result

def create_consolidated_vcf(imported, phased):
    """
    Create consolidated VCF.

    Parameters
    ----------
    imported : pypgx.Archive
        Archive file with the semantic type VcfFrame[Imported].
    phased : pypgx.Archive
        Archive file with the semandtic type VcfFrame[Phased].

    Returns
    -------
    pypgx.Archive
        Archive file with the semantic type VcfFrame[Consolidated].
    """
    vcf1 = sdk.Archive.from_file(imported)
    vcf2 = sdk.Archive.from_file(phased)

    vcf1.data = vcf1.data.strip('GT:AD:DP')
    vcf2.data = vcf2.data.strip('GT')

    def one_row(r):
        variant = f'{r.CHROM}-{r.POS}-{r.REF}-{r.ALT}'
        s = vcf1.data.fetch(variant)

        def one_gt(g):
            return ':'.join(g.split(':')[1:])

        s[9:] = s[9:].apply(one_gt)
        r[9:] = r[9:].str.cat(s[9:], sep=':')

        return r

    vf1 = pyvcf.VcfFrame([], vcf2.data.df.apply(one_row, axis=1))
    vf1.df.FORMAT = 'GT:AD:DP'

    vf2 = vcf1.data.filter_vcf(vcf2.data, opposite=True)
    vf3 = pyvcf.VcfFrame([], pd.concat([vf1.df, vf2.df])).sort()
    vf3 = vf3.pseudophase()

    metadata = vcf2.copy_metadata()
    metadata['SemanticType'] = 'VcfFrame[Consolidated]'

    result = sdk.Archive(metadata, vf3)

    return result

def estimate_phase_beagle(
    target, panel, impute=False
):
    """
    Estimate haplotype phase of observed variants with the Beagle program.

    Parameters
    ----------
    target : str
        Archive file with the semantic type VcfFrame[Imported].
    panel : str
        Reference haplotype panel.
    impute : bool, default: False
        Whether to perform imputation of missing genotypes.

    Returns
    -------
    pypgx.Archive
        Archive file with the semantic type VcfFrame[Phased].
    """
    vcf = sdk.Archive.from_file(target)
    region = get_region(vcf.metadata['Gene'], assembly=vcf.metadata['Assembly'])
    path = os.path.dirname(os.path.abspath(__file__))
    program = f'{path}/beagle.28Jun21.220.jar'
    metadata = vcf.copy_metadata()
    metadata['SemanticType'] = 'VcfFrame[Phased]'
    metadata['Program'] = 'Beagle'
    with tempfile.TemporaryDirectory() as t:
        vcf.data.to_file(f'{t}/input.vcf')
        command = [
            'java', '-Xmx2g', '-jar', program,
            f'gt={t}/input.vcf',
            f'chrom={region}',
            f'ref={panel}',
            f'out={t}/output',
            f'impute={str(impute).lower()}'
        ]
        subprocess.run(command)
        data = pyvcf.VcfFrame.from_file(f'{t}/output.vcf.gz')
    result = sdk.Archive(metadata, data)
    return result

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
    if gene not in list_genes():
        raise GeneNotFoundError(gene)

    if phenotype not in list_phenotypes():
        raise PhenotypeNotFoundError(phenotype)

    df = load_phenotype_table()
    i = (df.Gene == gene) & (df.Phenotype == phenotype)

    return df[i].Priority.values[0]

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

def import_vcf(gene, vcf, assembly='GRCh37'):
    """
    Import VCF data.
    """
    vf = pyvcf.VcfFrame.from_file(vcf)
    region = get_region(gene, assembly=assembly)
    data = vf.slice(region).strip('GT:AD:DP')
    metadata = {
        'Gene': gene,
        'Assembly': assembly,
        'SemanticType': 'VcfFrame[Imported]',
    }
    result = sdk.Archive(metadata, data)
    return result

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
    >>> pypgx.list_genes()[:5] # First five target genes
    ['CACNA1S', 'CFTR', 'CYP2A13', 'CYP2B6', 'CYP2C8']
    >>> pypgx.list_genes(mode='control')
    ['EGFR', 'RYR1', 'VDR']
    >>> pypgx.list_genes(mode='all')[:5] # First five genes in the table
    ['CACNA1S', 'CFTR', 'CYP2A13', 'CYP2B6', 'CYP2C8']
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

def load_control_table():
    """
    Load the control gene table.

    Returns
    -------
    pandas.DataFrame
        Requested table.

    Examples
    --------

    >>> import pypgx
    >>> df = pypgx.load_control_table()
    """
    b = BytesIO(pkgutil.get_data(__name__, 'data/control-table.csv'))
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

def predict_alleles(input):
    """
    Predict candidate star alleles based on observed variants.

    The input VCF must be fully phased.

    Parameters
    ----------
    input : pypgx.Archive
        Archive file with the semantic type VcfFrame[Consolidated].

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
    vcf = sdk.Archive.from_file(input)

    if vcf.metadata['Gene'] not in list_genes():
        raise GeneNotFoundError(vcf.metadata['Gene'])

    table = build_definition_table(vcf.metadata['Gene'],
        assembly=vcf.metadata['Assembly'])

    vf = vcf.data.filter_vcf(table)

    stars = {}

    func = lambda r: f'{r.CHROM}-{r.POS}-{r.REF}-{r.ALT}'

    for star in table.samples:
        df = table.df[table.df[star] == '1']
        stars[star] = set(df.apply(func, axis=1))

    samples = {}

    for sample in vf.samples:
        samples[sample] = [[], []]
        df = vf.df[sample].str.split(':').str[0].str.split('|', expand=True)
        df.index = vf.df.apply(func, axis=1)
        for i in [0, 1]:
            try:
                s = set(df[i][df[i] == '1'].index)
            except KeyError:
                s = set()
            for star, variants in stars.items():
                if variants.issubset(s):
                    samples[sample][i].append(star)
            if not samples[sample][i]:
                default = get_default_allele(vcf.metadata['Gene'],
                    vcf.metadata['Assembly'])
                if default:
                    samples[sample][i].append(default)
            samples[sample][i] = sort_alleles(vcf.metadata['Gene'], samples[sample][i])
            samples[sample][i] = ';'.join(samples[sample][i]) + ';'

    data = pd.DataFrame(samples).T
    data.columns = ['Haplotype1', 'Haplotype2']

    metadata = vcf.copy_metadata()
    metadata['SemanticType'] = 'TSV[Alleles]'
    result = sdk.Archive(metadata, data)
    return result

def predict_cnv(result):
    """
    Predict CNV based on copy number data.

    If there are missing values because, for example, the input data was
    generated with targeted sequencing, they will be filled in with the
    sample's median copy number.

    Parameters
    ----------
    result : pypgx.Archive or str
        Archive file with the semantic type CovFrame[CopyNumber].

    Returns
    -------
    pypgx.Archive
        Archive file with the semantic type TSV[CNVCalls].
    """
    result = sdk.Archive.from_file(result)
    path = os.path.dirname(os.path.abspath(__file__))
    model = sdk.Archive.from_file(f"{path}/cnv/{result.metadata['Gene']}.zip").data
    df = result.data.df.iloc[:, 2:]
    df = df.fillna(df.median())
    X = df.T.to_numpy()
    predictions = model.predict(X)
    df = load_cnv_table()
    cnvs = dict(zip(df.Code, df.Name))
    predictions = [cnvs[x] for x in predictions]
    metadata = result.copy_metadata()
    metadata['SemanticType'] = 'TSV[CNVCalls]'
    data = pd.DataFrame({'Sample': result.data.samples, 'CNV': predictions})
    return sdk.Archive(metadata, data)

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

def sort_alleles(gene, alleles, assembly='GRCh37'):
    """
    Sort star alleles by various creteria.

    Parameters
    ----------
    gene : str
        Gene name.
    alleles : str
        List of star allele.
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.

    Returns
    -------
    list
        Sorted list.

    Examples
    --------

    >>> import pypgx
    >>> pypgx.sort_alleles('CYP2D6', ['*1', '*4'])
    ['*4', '*1']
    """
    def f(allele):
        function = get_function(gene, allele)

        if pd.isna(function):
            a = len(FUNCTION_ORDER)
        else:
            a = FUNCTION_ORDER.index(function)

        variants = list_variants(gene, allele, assembly=assembly)

        b = len(variants) * -1

        return (a, b)

    return sorted(alleles, key=f)

def test_cnv_caller(caller, target, calls):
    """
    Test a CNV caller for target gene.

    Parameters
    ----------
    caller : pypgx.Archive
        Archive file with the semantic type Model[CNV].
    target : pypgx.Archive
        Archive file with the semantic type CovFrame[CopyNumber].
    calls : pypgx.Archive
        Archive file with the semantic type TSV[CNVCalls].
    """
    model = sdk.Archive.from_file(caller)
    copy_number = sdk.Archive.from_file(target)
    cnv_calls = sdk.Archive.from_file(calls)
    df = load_cnv_table()
    cnv_dict = dict(zip(df.Name, df.Code))
    cnv_calls.data['Code'] = cnv_calls.data.apply(lambda r: cnv_dict[r.CNV], axis=1)
    columns = ['Chromosome', 'Position'] + cnv_calls.data.Sample.to_list()
    copy_number.data.df = copy_number.data.df[columns]
    X = copy_number.data.df.iloc[:, 2:].T.to_numpy()
    Y = cnv_calls.data['Code'].to_numpy()
    predictions = model.data.predict(X)
    print(f'Accuracy: {sum(predictions == Y) / len(Y)} ({sum(predictions == Y)}/{len(Y)})')

def train_cnv_caller(target, calls):
    """
    Train a CNV caller for target gene.

    This method will return a SVM-based multiclass classifier that implements
    the one-vs-rest stategy.

    Parameters
    ----------
    target : pypgx.Archive
        Archive file with the semantic type CovFrame[CopyNumber].
    calls : pypgx.Archive
        Archive file with the semantic type TSV[CNVCalls].

    Returns
    -------
    pypgx.Archive
        Archive file with the semantic type Model[CNV].
    """
    copy_number = sdk.Archive.from_file(target)
    cnv_calls = sdk.Archive.from_file(calls)
    df = load_cnv_table()
    cnv_dict = dict(zip(df.Name, df.Code))
    cnv_calls.data['Code'] = cnv_calls.data.apply(lambda r: cnv_dict[r.CNV], axis=1)
    columns = ['Chromosome', 'Position'] + cnv_calls.data.Sample.to_list()
    copy_number.data.df = copy_number.data.df[columns]
    X = copy_number.data.df.iloc[:, 2:].T.to_numpy()
    Y = cnv_calls.data['Code'].to_numpy()
    model = OneVsRestClassifier(SVC(random_state=1)).fit(X, Y)
    metadata = copy_number.copy_metadata()
    metadata['SemanticType'] = 'Model[CNV]'
    return sdk.Archive(metadata, model)
