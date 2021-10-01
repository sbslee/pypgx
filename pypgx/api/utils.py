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

class NotTargetGeneError(Exception):
    """Raise if specified gene is not one of the target genes."""

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

    df = df.fillna(method='ffill')
    df = df.fillna(method='bfill')

    df.iloc[:, 2:] = df.iloc[:, 2:].apply(lambda c: median_filter(c, size=1000), axis=0)

    if df.isnull().values.any():
        raise ValueError('Missing values detected')

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
      CHROM       POS         ID REF ALT QUAL FILTER      INFO FORMAT *2 *3
    0    19  15990431  rs2108622   C   T    .      .  VI=V433M     GT  0  1
    1    19  16008388  rs3093105   A   C    .      .   VI=W12G     GT  1  0
    >>> vf = pypgx.build_definition_table('CYP4F2', assembly='GRCh38')
    >>> vf.df
      CHROM       POS         ID REF ALT QUAL FILTER      INFO FORMAT *2 *3
    0    19  15879621  rs2108622   C   T    .      .  VI=V433M     GT  0  1
    1    19  15897578  rs3093105   A   C    .      .   VI=W12G     GT  1  0
    """
    if gene not in list_genes():
        raise GeneNotFoundError(gene)

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
            v1 = set(list_variants(gene, a, assembly=assembly, mode='core'))
            v2 = set(list_variants(gene, b, assembly=assembly, mode='core'))
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
    genotypes : str or pypgx.Archive
        Archive file or object with the semantic type SampleTable[Genotypes].
    alleles : str or pypgx.Archive
        Archive file or object with the semantic type SampleTable[Alleles].
    cnv_calls : str or pypgx.Archive
        Archive file or object with the semantic type SampleTable[CNVCalls].

    Returns
    -------
    pypgx.Archive
        Archive object with the semantic type SampleTable[Results].
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
        Archive object with the semantic type SampleTable[Statistcs].
    """
    bam_files, bam_prefix = sdk.parse_input_bams(bam=bam, fn=fn)

    df = load_gene_table()

    if gene is not None:
        region = df[df.Gene == gene][f'{assembly}Region'].values[0]

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

def compute_copy_number(read_depth, control_statistics, samples=None):
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
    read_depth : str or pypgx.Archive
        Archive file or object with the semantic type CovFrame[ReadDepth].
    control_statistcs : str or pypgx.Archive
        Archive file or object with the semandtic type SampleTable[Statistics].
    samples : list, optional
        List of known samples without SV.

    Returns
    -------
    pypgx.Archive
        Archive file with the semandtic type CovFrame[CopyNumber].
    """
    if isinstance(read_depth, str):
        read_depth = sdk.Archive.from_file(read_depth)

    read_depth.check('CovFrame[ReadDepth]')

    if isinstance(control_statistics, str):
        control_statistics = sdk.Archive.from_file(control_statistics)

    control_statistics.check('SampleTable[Statistics]')

    if set(read_depth.data.samples) != set(control_statistics.data.index):
        raise ValueError('Different sample sets found')

    # Apply intra-sample normalization.
    df = read_depth.data.copy_df()
    medians = control_statistics.data['50%']
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
    metadata['Control'] = control_statistics.metadata['Control']
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
        Archive object with the semantic type CovFrame[ReadDepth].
    """
    metadata = {
        'Gene': gene,
        'Assembly': assembly,
        'SemanticType': 'CovFrame[ReadDepth]',
    }

    bam_files, bam_prefix = sdk.parse_input_bams(bam=bam, fn=fn)

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

    archive = sdk.Archive(metadata, data)

    return archive

def create_consolidated_vcf(imported_variants, phased_variants):
    """
    Create consolidated VCF.

    Parameters
    ----------
    imported_variants : str or pypgx.Archive
        Archive file or object with the semantic type VcfFrame[Imported].
    phased_variants : str or pypgx.Archive
        Archive file or object with the semandtic type VcfFrame[Phased].

    Returns
    -------
    pypgx.Archive
        Archive object with the semantic type VcfFrame[Consolidated].
    """
    if isinstance(imported_variants, str):
        imported_variants = sdk.Archive.from_file(imported_variants)

    imported_variants.check('VcfFrame[Imported]')

    if isinstance(phased_variants, str):
        phased_variants = sdk.Archive.from_file(phased_variants)

    phased_variants.check('VcfFrame[Phased]')

    if imported_variants.metadata['Gene'] != phased_variants.metadata['Gene']:
        raise ValueError('Found two different genes')

    gene = imported_variants.metadata['Gene']

    if imported_variants.metadata['Assembly'] != phased_variants.metadata['Assembly']:
        raise ValueError('Found two different assemblies')

    assembly = imported_variants.metadata['Assembly']

    vf1 = imported_variants.data.strip('GT:AD:DP:AF')
    vf2 = phased_variants.data.strip('GT')

    def one_row(r):
        variant = f'{r.CHROM}-{r.POS}-{r.REF}-{r.ALT}'
        s = vf1.fetch(variant)

        def one_gt(g):
            return ':'.join(g.split(':')[1:])

        s[9:] = s[9:].apply(one_gt)
        r[9:] = r[9:].str.cat(s[9:], sep=':')

        return r

    vf3 = pyvcf.VcfFrame([], vf2.df.apply(one_row, axis=1))
    vf3.df.INFO = 'Phased'
    vf3.df.FORMAT = 'GT:AD:DP:AF'

    vf4 = vf1.filter_vcf(vf2, opposite=True)
    vf5 = pyvcf.VcfFrame([], pd.concat([vf3.df, vf4.df])).sort()

    anchors = {}

    for i, r in vf2.df.iterrows():
        for allele in r.ALT.split(','):
            variant = f'{r.CHROM}-{r.POS}-{r.REF}-{allele}'
            for sample in vf2.samples:
                if sample not in anchors:
                    anchors[sample] = [[], []]
                gt = r[sample].split(':')[0].split('|')
                if gt[0] != '0':
                    anchors[sample][0].append(variant)
                if gt[1] != '0':
                    anchors[sample][1].append(variant)

    def one_row(r):
        if 'Phased' in r.INFO:
            return r

        r.FORMAT += ':PE'

        for sample in vf5.samples:
            if not pyvcf.gt_het(r[sample]):
                r[sample] = pyvcf.gt_pseudophase(r[sample]) + ':0,0,0,0'
                continue

            scores = [[0, 0], [0, 0]]

            gt = r[sample].split(':')[0].split('/')

            for i in [0, 1]:
                if gt[i] == '0':
                    continue

                alt_allele = r.ALT.split(',')[int(gt[i]) - 1]

                variant = f'{r.CHROM}-{r.POS}-{r.REF}-{alt_allele}'

                star_alleles = list_alleles(gene, variants=variant, assembly=assembly)

                for j in [0, 1]:
                    for star_allele in star_alleles:
                        score = 0
                        for x in anchors[sample][j]:
                            if x in list_variants(gene, star_allele, assembly=assembly, mode='all'):
                                score += 1
                        if score > scores[i][j]:
                            scores[i][j] = score

            a = scores[0][0]
            b = scores[0][1]
            c = scores[1][0]
            d = scores[1][1]

            if max([a, b]) == max([c, d]):
                if a < b and c > d:
                    flip = True
                elif a == b and c > d:
                    flip = True
                elif a < b and c == d:
                    flip = True
                else:
                    flip = False
            else:
                if max([a, b]) > max([c, d]):
                    if a > b:
                        flip = False
                    else:
                        flip = True
                else:
                    if c > d:
                        flip = True
                    else:
                        flip = False

            if flip:
                result = f'{gt[1]}|{gt[0]}'
            else:
                result = f'{gt[0]}|{gt[1]}'

            result = result + ':' + ':'.join(r[sample].split(':')[1:])
            r[sample] = result + ':' + ','.join([str(x) for x in scores[0] + scores[1]])

        return r

    vf5.df = vf5.df.apply(one_row, axis=1)

    metadata = phased_variants.copy_metadata()
    metadata['SemanticType'] = 'VcfFrame[Consolidated]'

    result = sdk.Archive(metadata, vf5)

    return result

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
    imported_variants, panel=None, impute=False
):
    """
    Estimate haplotype phase of observed variants with the Beagle program.

    Parameters
    ----------
    imported_variants : str or pypgx.Archive
        Archive file or object with the semantic type VcfFrame[Imported].
    panel : str, optional
        Reference haplotype panel. By default, the 1KGP panel is used.
    impute : bool, default: False
        Whether to perform imputation of missing genotypes.

    Returns
    -------
    pypgx.Archive
        Archive object with the semantic type VcfFrame[Phased].
    """
    if isinstance(imported_variants, str):
        imported_variants = sdk.Archive.from_file(imported_variants)

    imported_variants.check('VcfFrame[Imported]')

    gene = imported_variants.metadata['Gene']
    assembly = imported_variants.metadata['Assembly']
    region = get_region(gene, assembly=assembly)
    beagle = f'{PROGRAM_PATH}/pypgx/api/beagle.28Jun21.220.jar'
    if panel is None:
        panel = f'{PROGRAM_PATH}/pypgx/api/1kgp/{assembly}/{gene}.vcf.gz'

    metadata = imported_variants.copy_metadata()
    metadata['SemanticType'] = 'VcfFrame[Phased]'
    metadata['Program'] = 'Beagle'

    if imported_variants.data.empty:
        return sdk.Archive(metadata, imported_variants.data.copy())
    with tempfile.TemporaryDirectory() as t:
        imported_variants.data.to_file(f'{t}/input.vcf')
        command = [
            'java', '-Xmx2g', '-jar', beagle,
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
    archive : str or pypgx.archive
        Archive file or object.
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
    gene, depth_of_coverage
):
    """
    Import read depth data for target gene.

    Parameters
    ----------
    gene : str
        Gene name.
    depth_of_coverage : str or pypgx.Archive
        Archive file or object with the semantic type
        CovFrame[DepthOfCoverage].

    Returns
    -------
    pypgx.Archive
        Archive object with the semantic type CovFrame[ReadDepth].
    """
    if isinstance(depth_of_coverage, str):
        depth_of_coverage = sdk.Archive.from_file(depth_of_coverage)

    depth_of_coverage.check('CovFrame[DepthOfCoverage]')

    metadata = depth_of_coverage.copy_metadata()
    region = get_region(gene, assembly=metadata['Assembly'])
    data = depth_of_coverage.data.chr_prefix().slice(region)
    metadata['Gene'] = gene
    metadata['SemanticType'] = 'CovFrame[ReadDepth]'

    return sdk.Archive(metadata, data)

def import_variants(gene, vcf, assembly='GRCh37'):
    """
    Import variant data for the target gene.

    Parameters
    ----------
    gene : str
        Target gene.
    vcf : fuc.pyvcf.VcfFrame or str
        VCF file (zipped or unzipped).
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.

    Returns
    -------
    pypgx.Archive
        Archive object with the semantic type VcfFrame[Imported].
    """
    if isinstance(vcf, str):
        vf = pyvcf.VcfFrame.from_file(vcf)
    else:
        vf = vcf

    region = get_region(gene, assembly=assembly)

    data = vf.chr_prefix().slice(region).strip('GT:AD:DP')
    data = data.add_af().unphase()

    metadata = {
        'Gene': gene,
        'Assembly': assembly,
        'SemanticType': 'VcfFrame[Imported]',
    }

    return sdk.Archive(metadata, data)

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
    if gene not in list_genes():
        raise GeneNotFoundError(gene)

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

def list_variants(gene, allele, mode='all', assembly='GRCh37'):
    """
    List all variants that define specified allele.

    Some alleles will return an empty list because they do not contain any
    variants (e.g. reference allele).

    Parameters
    ----------
    gene : str
        Target gene.
    allele : str
        Star allele.
    mode : {'all', 'core', 'tag'}, default: 'all'
        Whether to return all variants, core variants only, or tag variants
        only.
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
    >>> pypgx.list_variants('CYP2B6', '*6', mode='all')
    ['19-41495755-T-C', '19-41496461-T-C', '19-41512841-G-T', '19-41515263-A-G']
    >>> pypgx.list_variants('CYP2B6', '*6', mode='core')
    ['19-41512841-G-T', '19-41515263-A-G']
    >>> pypgx.list_variants('CYP2B6', '*6', mode='tag')
    ['19-41495755-T-C', '19-41496461-T-C']
    """
    if gene not in list_genes():
        raise GeneNotFoundError(gene)

    df = load_allele_table()
    df = df[(df.Gene == gene) & (df.StarAllele == allele)]

    if df.empty:
        raise AlleleNotFoundError(gene, allele)

    core = df[f'{assembly}Core'].values[0]
    tag = df[f'{assembly}Tag'].values[0]

    if pd.isna(core):
        core = []
    else:
        core = core.split(',')

    if pd.isna(tag):
        tag = []
    else:
        tag = tag.split(',')

    if mode == 'all':
        results = tag + core
    elif mode == 'core':
        results = core
    elif mode == 'tag':
        results = tag
    else:
        raise ValueError('Incorrect mode')

    return results

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
                all_alleles += [x for x in candidates if x not in all_alleles]
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

def predict_cnv(copy_number, cnv_caller=None):
    """
    Predict CNV for target gene based on copy number data.

    If there are missing values because, for example, the input data was
    generated with targeted sequencing, they will be imputed with forward
    filling.

    Parameters
    ----------
    copy_number : str or pypgx.Archive
        Archive file or object with the semantic type CovFrame[CopyNumber].
    cnv_caller : str or pypgx.Archive, optional
        Archive file or object with the semantic type Model[CNV]. By default,
        a pre-trained CNV caller will be used.

    Returns
    -------
    pypgx.Archive
        Archive object with the semantic type SampleTable[CNVCalls].
    """
    if isinstance(copy_number, str):
        copy_number = sdk.Archive.from_file(copy_number)

    copy_number.check('CovFrame[CopyNumber]')

    if cnv_caller is None:
        cnv_caller = sdk.Archive.from_file(f"{PROGRAM_PATH}/pypgx/api/cnv/{copy_number.metadata['Gene']}.zip")
    else:
        if isinstance(cnv_caller, str):
            cnv_caller = sdk.Archive.from_file(cnv_caller)

        cnv_caller.check('Model[CNV]')

    copy_number = _process_copy_number(copy_number)

    df = copy_number.data.df.iloc[:, 2:]
    X = df.T.to_numpy()
    predictions = cnv_caller.data.predict(X)
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

def prepare_depth_of_coverage(
    bam=None, fn=None, assembly='GRCh37', bed=None
):
    """
    Prepare a depth of coverage file for all target genes with SV.

    By default, the input data is assumed to be WGS. If it's targeted
    sequencing, you must provide a BED file with ``bed`` to indicate
    probed regions.

    Parameters
    ----------
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
        Archive object with the semantic type CovFrame[DepthOfCoverage].
    """
    metadata = {
        'Assembly': assembly,
        'SemanticType': 'CovFrame[DepthOfCoverage]',
    }

    bam_files, bam_prefix = sdk.parse_input_bams(bam=bam, fn=fn)

    regions = create_regions_bed(
        merge=True, sv_genes=True
    ).gr.df.apply(
        lambda r: f'{r.Chromosome}:{r.Start}-{r.End}', axis=1
    ).to_list()

    cfs = []

    for region in regions:
        cf = pycov.CovFrame.from_bam(
            bam=bam_files, region=f'{bam_prefix}{region}', zero=True
        )
        cfs.append(cf)

    cf = pycov.concat(cfs)

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

    return sdk.Archive(metadata, cf)

def print_metadata(input):
    """
    Print the metadata of specified archive.

    Parameters
    ----------
    input : pypgx.Archive
        Archive file.
    """
    zf = zipfile.ZipFile(input)
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
        core = list_variants(gene, allele, assembly=assembly, mode='core')
        b = len(core) * -1
        return (a, b)

    return sorted(alleles, key=f)

def test_cnv_caller(
    cnv_caller, copy_number, cnv_calls, confusion_matrix=None
):
    """
    Test a CNV caller for the target gene.

    Parameters
    ----------
    cnv_caller : str or pypgx.Archive
        Archive file or object with the semantic type Model[CNV].
    copy_number : str or pypgx.Archive
        Archive file or object with the semantic type CovFrame[CopyNumber].
    cnv_calls : str or pypgx.Archive
        Archive file or object with the semantic type SampleTable[CNVCalls].
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
