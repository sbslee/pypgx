"""
The utils submodule is the main suite of tools for PGx research.
"""

import pkgutil
from io import BytesIO
import tempfile
import zipfile
import subprocess
import os
import pickle
import pathlib

from .. import sdk

import numpy as np
import pandas as pd
from fuc import pybam, pyvcf, pycov, common, pybed
from sklearn import model_selection, metrics
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import SVC
from sklearn.impute import SimpleImputer
from scipy.ndimage import median_filter

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

class PhenotypeNotFoundError(Exception):
    """Raise if specified phenotype is not present in the phenotype table."""

###################
# Private methods #
###################

def _process_copy_number(copy_number):
    df = copy_number.data.copy_df()
    region = get_region(copy_number.metadata['Gene'], assembly=copy_number.metadata['Assembly'])
    chrom, start, end = common.parse_region(region)
    if (end - start + 1) > copy_number.data.shape[0]:
        temp = pd.DataFrame.from_dict({'Temp': range(int(df.Position.iat[0]-1), int(df.Position.iat[-1])+1)})
        temp = temp.merge(df, left_on='Temp', right_on='Position', how='outer')
        df = temp.drop(columns='Temp')
    df = df.fillna(df.median())
    df.iloc[:, 2:] = df.iloc[:, 2:].apply(lambda c: median_filter(c, size=1000), axis=0)
    return sdk.Archive(copy_number.copy_metadata(), pycov.CovFrame(df))

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
            v1 = set(list_variants(gene, a, assembly=assembly))
            v2 = set(list_variants(gene, b, assembly=assembly))
            if v1.issubset(v2):
                result = True
                break
        results.append(result)

    return [x for i, x in enumerate(alleles) if not results[i]]

def combine_results(genotypes=None, alleles=None, cnv_calls=None):
    """
    Combine various results for the target gene.

    Parameters
    ----------
    genotypes : pypgx.Archive or str
        Archive file with the semantic type SampleTable[Genotypes].
    alleles : pypgx.Archive or str
        Archive file with the semantic type SampleTable[Alleles].
    cnv_calls : pypgx.Archive or str
        Archive file with the semantic type SampleTable[CNVCalls].

    Returns
    -------
    pypgx.Archive
        Archive file with the semantic type SampleTable[Results].
    """
    if isinstance(genotypes, str):
        genotypes = sdk.Archive.from_file(genotypes)

    if genotypes is not None:
        genotypes.check('SampleTable[Genotypes]')

    if isinstance(alleles, str):
        alleles = sdk.Archive.from_file(alleles)

    if alleles is not None:
        alleles.check('SampleTable[Alleles]')

    if isinstance(cnv_calls, str):
        cnv_calls = sdk.Archive.from_file(cnv_calls)

    if cnv_calls is not None:
        cnv_calls.check('SampleTable[CNVCalls]')

    tables = [x for x in [genotypes, alleles, cnv_calls] if x is not None]

    if not tables:
        raise ValueError('No input data detected')

    metadata = {}

    for k in ['Gene', 'Assembly']:
        l = [x.metadata[k] for x in tables]
        if len(set(l)) > 1:
            raise ValueError(f'Found incompatible inputs: {l}')
        metadata[k] = l[0]

    data = [x.data for x in tables]

    df = pd.concat(data, axis=1)

    cols = ['Genotype', 'Haplotype1', 'Haplotype2', 'AlternativePhase', 'VariantData', 'CNV']

    for col in cols:
        if col not in df.columns:
            df[col] = np.nan

    metadata['SemanticType'] = 'SampleTable[Results]'

    return sdk.Archive(metadata, df[cols])

def compute_control_statistics(
    bam=None, fn=None, gene=None, region=None, assembly='GRCh37', bed=None
):
    """
    Compute copy number from read depth for target gene.

    Input BAM files must be specified with either ``bam`` or ``fn``, but
    it's an error to use both. Similarly, control gene must be specified with
    either ``gene`` or ``region``, but it's an error to use both.

    By default, the input data is assumed to be WGS. If it's targeted
    sequencing, you must provide a BED file with ``bed`` to indicate
    probed regions.

    Parameters
    ----------
    bam : list, optional
        One or more BAM files.
    fn : str, optional
        File containing one BAM file per line.
    gene : str, optional
        Control gene (recommended choices: 'EGFR', 'RYR1', 'VDR').
    region : str, optional
        Custom region to use as control gene ('chrom:start-end').
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.
    bed : str, optional
        BED file.

    Returns
    -------
    pypgx.Archive
        Archive file with the semantic type SampleTable[Statistcs].
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
        bam_prefix = 'chr'
    else:
        bam_prefix = ''

    cf = pycov.CovFrame.from_bam(
        bam=bam_files, region=f'{bam_prefix}{region}', zero=False
    )

    metadata = {
        'Control': gene,
        'Assembly': assembly,
        'SemanticType': 'SampleTable[Statistics]',
    }

    if bed:
        metadata['Platform'] = 'Targeted'
        bf = pybed.BedFrame.from_file(bed)
        if any(['chr' in x for x in bf.contigs]):
            bed_prefix = 'chr'
        else:
            bed_prefix = ''
        if bam_prefix and bed_prefix:
            pass
        elif not bam_prefix and not bed_prefix:
            pass
        elif bam_prefix and not bed_prefix:
            bf = bf.chr_prefix(mode='add')
        else:
            bf = bf.chr_prefix(mode='remove')
        cf = cf.mask_bed(bf, opposite=True)
    else:
        metadata['Platform'] = 'WGS'

    data = cf.df.iloc[:, 2:].describe().T
    result = sdk.Archive(metadata, data)

    return result

def compute_copy_number(read_depth, control_statistcs, samples=None):
    """
    Compute copy number from read depth for target gene.

    The method will convert read depth from target gene to copy number by
    performing intra-sample normalization using summary statistics from
    control gene.

    If the input data was generated with targeted sequencing as opposed to
    WGS, the method will also apply inter-sample normalization using summary
    statistics across all samples. For best results, it is recommended to
    manually specify a list of known reference samples that do not have SV.

    Parameters
    ----------
    read_depth : pypgx.Archive
        Archive file with the semantic type CovFrame[ReadDepth].
    control_statistcs : pypgx.Archive
        Archive file with the semandtic type SampleTable[Statistcs].
    samples : list, optional
        List of known samples without SV.

    Returns
    -------
    pypgx.Archive
        Archive file with the semandtic type CovFrame[CopyNumber].
    """
    if isinstance(read_depth, str):
        read_depth = sdk.Archive.from_file(read_depth)

    if isinstance(control_statistcs, str):
        control_statistcs = sdk.Archive.from_file(control_statistcs)

    # Apply intra-sample normalization.
    df = read_depth.data.copy_df()
    medians = control_statistcs.data['50%']
    df.iloc[:, 2:] = df.iloc[:, 2:] / medians * 2

    # Apply inter-sample normalization.
    if read_depth.metadata['Platform'] == 'Targeted':
        if samples is None:
            medians = df.iloc[:, 2:].median(axis=1).replace(0, np.nan)
        else:
            medians = df[samples].median(axis=1).replace(0, np.nan)
        df.iloc[:, 2:] = df.iloc[:, 2:].div(medians, axis=0) * 2

    cf = pycov.CovFrame(df)
    metadata = read_depth.copy_metadata()
    metadata['SemanticType'] = 'CovFrame[CopyNumber]'
    metadata['Control'] = control_statistcs.metadata['Control']
    if samples is None:
        metadata['Samples'] = 'None'
    else:
        metadata['Samples'] = ','.join(samples)

    return sdk.Archive(metadata, cf)

def compute_target_depth(
    gene, bam=None, fn=None, assembly='GRCh37', bed=None
):
    """
    Compute read depth for target gene with BAM data.

    Input BAM files must be specified with either ``bam`` or ``fn``, but
    it's an error to use both.

    By default, the input data is assumed to be WGS. If it's targeted
    sequencing, you must provide a BED file with ``bed`` to indicate
    probed regions.

    Parameters
    ----------
    gene : str
        Target gene.
    bam : list, optional
        One or more BAM files.
    fn : str, optional
        File containing one BAM file per line.
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.
    bed : str, optional
        BED file.

    Returns
    -------
    pypgx.Archive
        Archive file with the semantic type CovFrame[ReadDepth].
    """
    metadata = {
        'Gene': gene,
        'Assembly': assembly,
        'SemanticType': 'CovFrame[ReadDepth]',
    }

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
        bam_prefix = 'chr'
    else:
        bam_prefix = ''

    region = get_region(gene, assembly=assembly)

    data = pycov.CovFrame.from_bam(
        bam=bam_files, region=f'{bam_prefix}{region}', zero=True
    )

    if bed:
        metadata['Platform'] = 'Targeted'
        bf = pybed.BedFrame.from_file(bed)
        if any(['chr' in x for x in bf.contigs]):
            bed_prefix = 'chr'
        else:
            bed_prefix = ''
        if bam_prefix and bed_prefix:
            pass
        elif not bam_prefix and not bed_prefix:
            pass
        elif bam_prefix and not bed_prefix:
            bf = bf.chr_prefix(mode='add')
        else:
            bf = bf.chr_prefix(mode='remove')
        data = data.mask_bed(bf, opposite=True)
    else:
        metadata['Platform'] = 'WGS'

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
    if isinstance(imported, str):
        imported = sdk.Archive.from_file(imported)

    if isinstance(phased, str):
        phased = sdk.Archive.from_file(phased)

    imported.data = imported.data.strip('GT:AD:DP:AF')
    phased.data = phased.data.strip('GT')

    def one_row(r):
        variant = f'{r.CHROM}-{r.POS}-{r.REF}-{r.ALT}'
        s = imported.data.fetch(variant)

        def one_gt(g):
            return ':'.join(g.split(':')[1:])

        s[9:] = s[9:].apply(one_gt)
        r[9:] = r[9:].str.cat(s[9:], sep=':')

        return r

    vf1 = pyvcf.VcfFrame([], phased.data.df.apply(one_row, axis=1))
    vf1.df.FORMAT = 'GT:AD:DP:AF'

    vf2 = imported.data.filter_vcf(phased.data, opposite=True)
    vf3 = pyvcf.VcfFrame([], pd.concat([vf1.df, vf2.df])).sort()
    vf3 = vf3.pseudophase()

    metadata = phased.copy_metadata()
    metadata['SemanticType'] = 'VcfFrame[Consolidated]'

    result = sdk.Archive(metadata, vf3)

    return result

def create_read_depth_tsv(bam=None, fn=None, assembly='GRCh37'):
    """
    Create TSV file containing read depth for genes with SV.

    Parameters
    ----------
    bam : list, optional
        One or more BAM files.
    fn : str, optional
        File containing one BAM file per line.
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.

    Returns
    -------
    fuc.pycov.CovFrame
        CovFrame object.
    """
    regions = create_regions_bed(
        merge=True, sv_genes=True
    ).gr.df.apply(
        lambda r: f'{r.Chromosome}:{r.Start}-{r.End}', axis=1
    ).to_list()

    cfs = []

    for region in regions:
        cf = pycov.CovFrame.from_bam(bam=bam, fn=fn, region=region, zero=True)
        cfs.append(cf)

    return pycov.concat(cfs)

def create_regions_bed(
    assembly='GRCh37', chr_prefix=False, merge=False, sv_genes=False
):
    """
    Create a BED file which contains all regions used by PyPGx.

    Parameters
    ----------
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.
    chr_prefix : bool, default: False
        Whether to add the 'chr' string in contig names.
    merge : bool, default: False
        Whether to merge overlapping intervals (gene names will be removed
        too).
    sv_genes : bool, default: False
        Whether to only return genes with SV.

    Returns
    -------
    fuc.pybed.BedFrame
        BED file.
    """
    df = load_gene_table()
    if sv_genes:
        df = df[df.SV]
    data = []
    for i, r in df.iterrows():
        region = r[f'{assembly}Region']
        fields = list(common.parse_region(region))
        fields.append(r.Gene)
        data.append(fields)
    df = pd.DataFrame(data)
    df.columns = ['Chromosome', 'Start', 'End', 'Name']
    bf = pybed.BedFrame.from_frame([], df)
    if chr_prefix:
        bf = bf.chr_prefix(mode='add')
    if merge:
        bf = bf.merge()
    return bf

def estimate_phase_beagle(
    imported_variants, panel, impute=False
):
    """
    Estimate haplotype phase of observed variants with the Beagle program.

    If your input data is GRCh37, I recommend using the 1000 Genomes Project
    phase 3 reference panel. You can easily download it thanks to the authors
    of Beagle:

    .. code-block:: console

        $ wget -r --no-parent http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/

    Parameters
    ----------
    imported_variants : str
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
    if isinstance(imported_variants, str):
        imported_variants = sdk.Archive.from_file(imported_variants)

    imported_variants.check('VcfFrame[Imported]')

    region = get_region(imported_variants.metadata['Gene'], assembly=imported_variants.metadata['Assembly'])
    path = os.path.dirname(os.path.abspath(__file__))
    program = f'{path}/beagle.28Jun21.220.jar'
    metadata = imported_variants.copy_metadata()
    metadata['SemanticType'] = 'VcfFrame[Phased]'
    metadata['Program'] = 'Beagle'
    if imported_variants.data.empty:
        return sdk.Archive(metadata, imported_variants.data.copy())
    with tempfile.TemporaryDirectory() as t:
        imported_variants.data.to_file(f'{t}/input.vcf')
        command = [
            'java', '-Xmx2g', '-jar', program,
            f'gt={t}/input.vcf',
            f'chrom={region}',
            f'ref={panel}',
            f'out={t}/output',
            f'impute={str(impute).lower()}'
        ]
        subprocess.run(command, stdout=subprocess.DEVNULL)
        data = pyvcf.VcfFrame.from_file(f'{t}/output.vcf.gz')
    return sdk.Archive(metadata, data)

def filter_samples(archive, samples=None, exclude=False, fn=None):
    """
    Filter Archive for specified samples.

    Samples can be specified with either ``samples`` or ``fn``, but it's an
    error to use both.

    Parameters
    ----------
    archive : pypgx.archive or str
        Archive object or file.
    samples : str or list
        Sample name or list of names (the order matters).
    exclude : bool, default: False
        If True, exclude specified samples.
    fn : str
        File containing one filename per line.

    Returns
    -------
    pypgx.Archive
        Fitlered Archive object.
    """
    if isinstance(archive, str):
        archive = sdk.Archive.from_file(archive)

    if isinstance(samples, str):
        samples = [samples]

    if samples is not None and fn is None:
        pass
    elif samples is None and fn is not None:
        samples = common.convert_file2list(fn)
    elif samples is not None and fn is not None:
        raise ValueError('Found two sets of samples')
    else:
        raise ValueError('Samples not found')

    if 'CovFrame' in archive.metadata['SemanticType']:
        data = archive.data.subset(samples, exclude=exclude)
    elif 'SampleTable' in archive.metadata['SemanticType']:
        data = archive.data.loc[samples]
    else:
        pass

    return sdk.Archive(archive.copy_metadata(), data)

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
    if gene not in list_genes():
        raise GeneNotFoundError(gene)

    if not has_score(gene):
        return np.nan

    df = load_allele_table()
    df = df[(df.Gene == gene) & (df.StarAllele == allele)]

    if df.empty:
        raise AlleleNotFoundError(gene + allele)

    return df.ActivityScore.values[0]

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

def import_read_depth(
    gene, read_depth, assembly='GRCh37', platform='WGS'
):
    """
    Import read depth data for target gene.

    Parameters
    ----------
    gene : str
        Gene name.
    read_depth : fuc.pycov.CovFrame or str
        TSV file containing read depth (zipped or unzipped).
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.
    platform : {'WGS', 'Targeted'}, default: 'WGS'
        NGS platform.

    Returns
    -------
    pypgx.Archive
        Archive file with the semantic type CovFrame[ReadDepth].
    """
    if isinstance(read_depth, str):
        cf = pycov.CovFrame.from_file(read_depth)
    else:
        cf = read_depth

    region = get_region(gene, assembly=assembly)

    data = cf.chr_prefix().slice(region)

    metadata = {
        'Gene': gene,
        'Assembly': assembly,
        'SemanticType': 'CovFrame[ReadDepth]',
        'Platform': platform,
    }

    return sdk.Archive(metadata, data)

def import_variants(gene, vcf, assembly='GRCh37'):
    """
    Import variant data for target gene.

    Parameters
    ----------
    gene : str
        Gene name.
    vcf : fuc.pyvcf.VcfFrame or str
        VCF file (zipped or unzipped).
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.

    Returns
    -------
    pypgx.Archive
        Archive file with the semantic type VcfFrame[Imported].
    """
    if isinstance(vcf, str):
        vf = pyvcf.VcfFrame.from_file(vcf)
    else:
        vf = vcf

    region = get_region(gene, assembly=assembly)

    data = vf.slice(region).strip('GT:AD:DP')
    data = data.add_af()

    metadata = {
        'Gene': gene,
        'Assembly': assembly,
        'SemanticType': 'VcfFrame[Imported]',
    }

    return sdk.Archive(metadata, data)

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

def predict_alleles(consolidated_variants):
    """
    Predict candidate star alleles based on observed SNVs and INDELs.

    Parameters
    ----------
    consolidated_variants : str or pypgx.Archive
        Archive file or object with the semantic type VcfFrame[Consolidated].

    Returns
    -------
    pypgx.Archive
        Archive object with the semantic type VcfFrame SampleTable[Alleles].

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
    if isinstance(consolidated_variants, str):
        consolidated_variants = sdk.Archive.from_file(consolidated_variants)
    consolidated_variants.check('VcfFrame[Consolidated]')
    gene = consolidated_variants.metadata['Gene']
    assembly = consolidated_variants.metadata['Assembly']
    definition_table = build_definition_table(gene, assembly)
    vf = consolidated_variants.data.filter_vcf(definition_table)
    ref_allele = get_ref_allele(gene, assembly)
    default_allele = get_default_allele(gene, assembly)

    stars = {}

    func = lambda r: f'{r.CHROM}-{r.POS}-{r.REF}-{r.ALT}'

    for star in definition_table.samples:
        df = definition_table.df[definition_table.df[star] == '1']
        stars[star] = set(df.apply(func, axis=1))

    samples = {}

    def one_haplotype(observed):
        candidates = []
        for star, variants in stars.items():
            if variants.issubset(observed):
                candidates.append(star)
        candidates = collapse_alleles(gene, candidates, assembly=assembly)
        if ref_allele != default_allele and ref_allele not in candidates and default_allele not in candidates:
            candidates.append(default_allele)
        if not candidates:
            candidates.append(default_allele)
        candidates = sort_alleles(gene, candidates)
        return candidates

    for sample in vf.samples:
        samples[sample] = []
        df = vf.df[sample].str.split(':').str[0].str.split('|', expand=True)
        df.index = vf.df.apply(func, axis=1)
        alt_phase = []
        all_alleles = []
        for i in [0, 1, 2]:
            try:
                observed = set(df[i][df[i] == '1'].index)
            except KeyError:
                observed = set()
            if i != 2:
                alt_phase += [x for x in observed if x not in alt_phase]
                candidates = one_haplotype(observed)
                all_alleles += [x for x in candidates if x not in all_alleles]
            else:
                candidates = one_haplotype(set(alt_phase))
                candidates = [x for x in candidates if x not in all_alleles]
                all_alleles = sort_alleles(gene, all_alleles)
            samples[sample].append(';'.join(candidates) + ';')

        af_list = []

        for allele in all_alleles:
            if allele == default_allele:
                af_list.append(f'{allele}:default')
            else:
                variants = ','.join(stars[allele])
                fractions = ','.join([str(vf.get_af(sample, x)) for x in stars[allele]])
                af_list.append(f'{allele}:{variants}:{fractions}')
        samples[sample].append(';'.join(af_list) + ';')

    data = pd.DataFrame(samples).T
    data.columns = ['Haplotype1', 'Haplotype2', 'AlternativePhase', 'VariantData']
    metadata = consolidated_variants.copy_metadata()
    metadata['SemanticType'] = 'SampleTable[Alleles]'
    result = sdk.Archive(metadata, data)
    return result

def predict_cnv(copy_number):
    """
    Predict CNV for target gene based on copy number data.

    If there are missing values because, for example, the input data was
    generated with targeted sequencing, they will be filled in with the
    sample's median copy number.

    Parameters
    ----------
    target : pypgx.Archive or str
        Archive file with the semantic type CovFrame[CopyNumber].

    Returns
    -------
    pypgx.Archive
        Archive file with the semantic type SampleTable[CNVCalls].
    """
    if isinstance(copy_number, str):
        copy_number = sdk.Archive.from_file(copy_number)
    copy_number = _process_copy_number(copy_number)
    path = os.path.dirname(os.path.abspath(__file__))
    model = sdk.Archive.from_file(f"{path}/cnv/{copy_number.metadata['Gene']}.zip").data
    df = copy_number.data.df.iloc[:, 2:]
    df = df.fillna(df.median())
    X = df.T.to_numpy()
    predictions = model.predict(X)
    df = load_cnv_table()
    df = df[df.Gene == copy_number.metadata['Gene']]
    cnvs = dict(zip(df.Code, df.Name))
    predictions = [cnvs[x] for x in predictions]
    metadata = copy_number.copy_metadata()
    metadata['SemanticType'] = 'SampleTable[CNVCalls]'
    data = pd.DataFrame({'CNV': predictions})
    data.index = copy_number.data.samples
    return sdk.Archive(metadata, data)

def predict_phenotype(gene, a, b):
    """
    Predict phenotype based on two haplotype calls.

    The method can handle star alleles with SV including gene deletion,
    duplication, and tandem arrangement.

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

def print_metadata(archive):
    """
    Print the metadata of specified archive.

    Parameters
    ----------
    archive : pypgx.Archive
        Archive file.
    """
    zf = zipfile.ZipFile(archive)
    parent = zf.filelist[0].filename.split('/')[0]
    with zf.open(f'{parent}/metadata.txt') as f:
        print(f.read().decode('utf-8').strip())

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
        a = FUNCTION_ORDER.index(function)
        variants = list_variants(gene, allele, assembly=assembly)
        b = len(variants) * -1
        return (a, b)

    return sorted(alleles, key=f)

def test_cnv_caller(
    cnv_caller, copy_number, cnv_calls, confusion_matrix=None
):
    """
    Test a CNV caller for the target gene.

    Parameters
    ----------
    cnv_caller : pypgx.Archive
        Archive file with the semantic type Model[CNV].
    copy_number : pypgx.Archive
        Archive file with the semantic type CovFrame[CopyNumber].
    cnv_calls : pypgx.Archive
        Archive file with the semantic type SampleTable[CNVCalls].
    confusion_matrix : str, optional
        Write the confusion matrix as a CSV file.
    """
    if isinstance(cnv_caller, str):
        cnv_caller = sdk.Archive.from_file(cnv_caller)

    cnv_caller.check('Model[CNV]')

    if isinstance(copy_number, str):
        copy_number = sdk.Archive.from_file(copy_number)

    copy_number.check('CovFrame[CopyNumber]')

    if isinstance(cnv_calls, str):
        cnv_calls = sdk.Archive.from_file(cnv_calls)

    cnv_calls.check('SampleTable[CNVCalls]')

    if not cnv_caller.metadata['Gene'] == copy_number.metadata['Gene'] == cnv_calls.metadata['Gene']:
        raise ValueError(f"Model[CNV] has {cnv_caller.metadata['Gene']}, CovFrame[CopyNumber] has {copy_number.metadata['Gene']}, and SampleTable[CNVCalls] has {cnv_calls.metadata['Gene']}")

    copy_number = _process_copy_number(copy_number)

    cnv_table = load_cnv_table()
    cnv_table = cnv_table[cnv_table.Gene == copy_number.metadata['Gene']]
    name2code = dict(zip(cnv_table.Name, cnv_table.Code))
    code2name = dict(zip(cnv_table.Code, cnv_table.Name))

    cnv_calls.data['Code'] = cnv_calls.data.apply(lambda r: name2code[r.CNV], axis=1)
    columns = ['Chromosome', 'Position'] + cnv_calls.data.index.to_list()
    copy_number.data.df = copy_number.data.df[columns]
    X = copy_number.data.df.iloc[:, 2:].T.to_numpy()
    Y = cnv_calls.data['Code'].to_numpy()
    predictions = cnv_caller.data.predict(X)
    results = predictions == Y
    print(f'Accuracy: {sum(results)/len(Y):.3f} ({sum(results)}/{len(Y)})')

    if confusion_matrix is not None:
        Y = [code2name[x] for x in Y]
        predictions = [code2name[x] for x in predictions]
        labels = cnv_table.Name.to_list()
        df = pd.DataFrame(metrics.confusion_matrix(Y, predictions, labels=labels))
        df.columns = labels
        df.index = labels
        df.to_csv(confusion_matrix)

def train_cnv_caller(copy_number, cnv_calls, confusion_matrix=None):
    """
    Train a CNV caller for the target gene.

    This method will return a SVM-based multiclass classifier that makes CNV
    calls using the one-vs-rest stategy.

    Parameters
    ----------
    copy_number : str or pypgx.Archive
        Archive file or object with the semantic type CovFrame[CopyNumber].
    cnv_calls : str or pypgx.Archive
        Archive file or object with the semantic type SampleTable[CNVCalls].
    confusion_matrix : str, optional
        Write the confusion matrix as a CSV file.

    Returns
    -------
    pypgx.Archive
        Archive object with the semantic type Model[CNV].
    """
    if isinstance(copy_number, str):
        copy_number = sdk.Archive.from_file(copy_number)

    copy_number.check('CovFrame[CopyNumber]')

    if isinstance(cnv_calls, str):
        cnv_calls = sdk.Archive.from_file(cnv_calls)

    cnv_calls.check('SampleTable[CNVCalls]')

    copy_number = _process_copy_number(copy_number)

    if copy_number.metadata['Gene'] != cnv_calls.metadata['Gene']:
        raise ValueError(f"CovFrame[CopyNumber] has {copy_number.metadata['Gene']}, while SampleTable[CNVCalls] has {cnv_calls.metadata['Gene']}")

    cnv_table = load_cnv_table()
    cnv_table = cnv_table[cnv_table.Gene == copy_number.metadata['Gene']]
    name2code = dict(zip(cnv_table.Name, cnv_table.Code))
    code2name = dict(zip(cnv_table.Code, cnv_table.Name))
    cnv_calls.data['Code'] = cnv_calls.data.apply(lambda r: name2code[r.CNV], axis=1)
    columns = ['Chromosome', 'Position'] + cnv_calls.data.index.to_list()
    copy_number.data.df = copy_number.data.df[columns]
    X = copy_number.data.df.iloc[:, 2:].T.to_numpy()
    Y = cnv_calls.data['Code'].to_numpy()
    model = OneVsRestClassifier(SVC(random_state=1)).fit(X, Y)
    metadata = copy_number.copy_metadata()
    metadata['SemanticType'] = 'Model[CNV]'
    predictions = model.predict(X)
    results = predictions == Y
    print(f'Accuracy: {sum(results)/len(Y):.3f} ({sum(results)}/{len(Y)})')

    if confusion_matrix is not None:
        Y = [code2name[x] for x in Y]
        predictions = [code2name[x] for x in predictions]
        labels = cnv_table.Name.to_list()
        df = pd.DataFrame(metrics.confusion_matrix(Y, predictions, labels=labels))
        df.columns = labels
        df.index = labels
        df.to_csv(confusion_matrix)

    return sdk.Archive(metadata, model)
