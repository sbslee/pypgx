"""
The utils submodule contains main actions of PyPGx.
"""

import pkgutil
from io import BytesIO
import tempfile
import zipfile
import subprocess
import os
import pickle

from . import core
from .. import sdk

import numpy as np
import pandas as pd
from fuc import pybam, pyvcf, pycov, common, pybed
from sklearn import model_selection, metrics
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import SVC
from sklearn.impute import SimpleImputer
from scipy.ndimage import median_filter

###################
# Private methods #
###################

def _process_copy_number(copy_number):
    df = copy_number.data.copy_df()
    region = core.get_region(copy_number.metadata['Gene'], assembly=copy_number.metadata['Assembly'])
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

def call_phenotypes(genotypes):
    """
    Call phenotypes for the target gene.

    Parameters
    ----------
    genotypes : str or pypgx.Archive
        Archive file or object with the semantic type SampleTable[Genotypes].

    Returns
    -------
    pypgx.Archive
        Archive object with the semantic type SampleTable[Phenotypes].
    """
    if isinstance(genotypes, str):
        genotypes = sdk.Archive.from_file(genotypes)

    genotypes.check_type('SampleTable[Genotypes]')

    gene = genotypes.metadata['Gene']

    def one_row(r):
        if r.Genotype == 'Indeterminate':
            phenotype = 'Indeterminate'
        else:
            a1, a2 = r.Genotype.split('/')
            phenotype = core.predict_phenotype(gene, a1, a2)
        return phenotype

    data = genotypes.data.apply(one_row, axis=1).to_frame()
    data.columns = ['Phenotype']

    metadata = {}
    metadata['Gene'] = gene
    metadata['SemanticType'] = 'SampleTable[Phenotypes]'

    return sdk.utils.Archive(metadata, data)

def combine_results(
    genotypes=None, phenotypes=None, alleles=None, cnv_calls=None
):
    """
    Combine various results for the target gene.

    Parameters
    ----------
    genotypes : str or pypgx.Archive, optional
        Archive file or object with the semantic type SampleTable[Genotypes].
    phenotypes : str or pypgx.Archive, optional
        Archive file or object with the semantic type SampleTable[Phenotypes].
    alleles : str or pypgx.Archive, optional
        Archive file or object with the semantic type SampleTable[Alleles].
    cnv_calls : str or pypgx.Archive, optional
        Archive file or object with the semantic type SampleTable[CNVCalls].

    Returns
    -------
    pypgx.Archive
        Archive object with the semantic type SampleTable[Results].
    """
    if isinstance(genotypes, str):
        genotypes = sdk.Archive.from_file(genotypes)

    if genotypes is not None:
        genotypes.check_type('SampleTable[Genotypes]')

    if isinstance(phenotypes, str):
        phenotypes = sdk.Archive.from_file(phenotypes)

    if phenotypes is not None:
        phenotypes.check_type('SampleTable[Phenotypes]')

    if isinstance(alleles, str):
        alleles = sdk.Archive.from_file(alleles)

    if alleles is not None:
        alleles.check_type('SampleTable[Alleles]')

    if isinstance(cnv_calls, str):
        cnv_calls = sdk.Archive.from_file(cnv_calls)

    if cnv_calls is not None:
        cnv_calls.check_type('SampleTable[CNVCalls]')

    tables = [x for x in [genotypes, phenotypes, alleles, cnv_calls]
        if x is not None]

    if not tables:
        raise ValueError('No input data detected')

    metadata = {}

    for k in ['Gene', 'Assembly']:
        l = [x.metadata[k] for x in tables if k in x.metadata]
        if len(set(l)) > 1:
            raise ValueError(f'Found incompatible inputs: {l}')
        metadata[k] = l[0]

    data = [x.data for x in tables]

    df = pd.concat(data, axis=1)

    cols = ['Genotype', 'Phenotype', 'Haplotype1', 'Haplotype2', 'AlternativePhase', 'VariantData', 'CNV']

    for col in cols:
        if col not in df.columns:
            df[col] = np.nan

    metadata['SemanticType'] = 'SampleTable[Results]'

    return sdk.Archive(metadata, df[cols])

def compare_genotypes(first, second, verbose=False):
    """
    Calculate concordance rate between two genotype results.

    The method will only use samples that appear in both genotype results.

    Parameters
    ----------
    first : str or pypgx.Archive
        First archive file or object with the semantic type
        SampleTable[Results].
    second : str or pypgx.Archive
        Second archive file or object with the semantic type
        SampleTable[Results].
    verbose : bool, default: False
        If True, print the verbose version of output.

    Examples
    --------
    >>> import pypgx
    >>> pypgx.compare_genotypes('first-results.zip', 'second-results.zip')
    Total: 32
    Compared: 32
    Concordance: 1.000 (32/32)
    """
    if isinstance(first, str):
        first = sdk.Archive.from_file(first)

    first.check_type('SampleTable[Results]')

    if isinstance(second, str):
        second = sdk.Archive.from_file(second)

    second.check_type('SampleTable[Results]')

    df = pd.concat([first.data.Genotype, second.data.Genotype], axis=1)

    print(f'Total: {df.shape[0]}')

    df.columns = ['First', 'Second']
    df = df.dropna()

    print(f'Compared: {df.shape[0]}')

    df['Concordant'] = df.First == df.Second

    print(f'Concordance: {sum(df.Concordant)/len(df.Concordant):.3f} ({sum(df.Concordant)}/{len(df.Concordant)})')

    if verbose:
        print('Discordant genotypes:')
        print(df[~df.Concordant])

def compute_control_statistics(
    bam=None, fn=None, gene=None, region=None, assembly='GRCh37', bed=None
):
    """
    Compute summary statistics for the control gene from BAM files.

    Parameters
    ----------
    bam : list, optional
        One or more BAM files. Cannot be used with ``fn``.
    fn : str, optional
        File containing one BAM file per line. Cannot be used with ``bam``.
    gene : str, optional
        Control gene (recommended choices: 'EGFR', 'RYR1', 'VDR'). Cannot be
        used with ``region``.
    region : str, optional
        Custom region to use as control gene ('chrom:start-end'). Cannot be
        used with ``gene``.
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.
    bed : str, optional
        By default, the input data is assumed to be WGS. If it targeted
        sequencing, you must provide a BED file to indicate probed regions.

    Returns
    -------
    pypgx.Archive
        Archive object with the semantic type SampleTable[Statistcs].
    """
    bam_files, bam_prefix = sdk.parse_input_bams(bam=bam, fn=fn)

    df = core.load_gene_table()

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
        if bf.has_chr_prefix:
            bed_prefix = 'chr'
        else:
            bed_prefix = ''
        if bam_prefix and bed_prefix:
            pass
        elif not bam_prefix and not bed_prefix:
            pass
        elif bam_prefix and not bed_prefix:
            bf = bf.update_chr_prefix(mode='add')
        else:
            bf = bf.update_chr_prefix(mode='remove')
        cf = cf.mask_bed(bf, opposite=True)
    else:
        metadata['Platform'] = 'WGS'

    data = cf.df.iloc[:, 2:].describe().T
    result = sdk.Archive(metadata, data)

    return result

def compute_copy_number(
    read_depth, control_statistics, samples_without_sv=None
):
    """
    Compute copy number from read depth for the target gene.

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
        Archive file or object with the semandtic type
        SampleTable[Statistics].
    samples_without_sv : list, optional
        List of known samples without SV.

    Returns
    -------
    pypgx.Archive
        Archive file with the semandtic type CovFrame[CopyNumber].
    """
    if isinstance(read_depth, str):
        read_depth = sdk.Archive.from_file(read_depth)

    read_depth.check_type('CovFrame[ReadDepth]')

    if isinstance(control_statistics, str):
        control_statistics = sdk.Archive.from_file(control_statistics)

    control_statistics.check_type('SampleTable[Statistics]')

    if set(read_depth.data.samples) != set(control_statistics.data.index):
        raise ValueError('Different sample sets found')

    # Apply intra-sample normalization.
    df = read_depth.data.copy_df()
    medians = control_statistics.data['50%']
    df.iloc[:, 2:] = df.iloc[:, 2:] / medians * 2

    # Apply inter-sample normalization.
    if read_depth.metadata['Platform'] == 'Targeted':
        if samples_without_sv is None:
            medians = df.iloc[:, 2:].median(axis=1).replace(0, np.nan)
        else:
            medians = df[samples_without_sv].median(axis=1).replace(0, np.nan)
        df.iloc[:, 2:] = df.iloc[:, 2:].div(medians, axis=0) * 2

    cf = pycov.CovFrame(df)
    metadata = read_depth.copy_metadata()
    metadata['SemanticType'] = 'CovFrame[CopyNumber]'
    metadata['Control'] = control_statistics.metadata['Control']
    if samples_without_sv is None:
        metadata['Samples'] = 'None'
    else:
        metadata['Samples'] = ','.join(samples_without_sv)

    return sdk.Archive(metadata, cf)

def compute_target_depth(
    gene, bam=None, fn=None, assembly='GRCh37', bed=None
):
    """
    Compute read depth for the target gene from BAM files.

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

    region = core.get_region(gene, assembly=assembly)

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
            bf = bf.update_chr_prefix(mode='add')
        else:
            bf = bf.update_chr_prefix(mode='remove')
        data = data.mask_bed(bf, opposite=True)
    else:
        metadata['Platform'] = 'WGS'

    archive = sdk.Archive(metadata, data)

    return archive

def count_alleles(results):
    """
    Count star alleles from genotype calls.
    """
    if isinstance(results, str):
        results = sdk.Archive.from_file(results)

    results.check_type('SampleTable[Results]')

    df = results.data.copy()

    def one_row(r):
        if r.Genotype == 'Indeterminate':
            return ['Indeterminate', 'Indeterminate']
        return r.Genotype.split('/')
    df = df.apply(one_row, axis=1, result_type='expand')
    s = pd.concat([df[0], df[1]])
    s = s.reset_index(drop=True)
    s = s.value_counts()
    s = s[core.sort_alleles(s.index.to_list(), by='name')]
    return s

def create_consolidated_vcf(imported_variants, phased_variants):
    """
    Create a consolidated VCF file.

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

    imported_variants.check_type('VcfFrame[Imported]')

    if isinstance(phased_variants, str):
        phased_variants = sdk.Archive.from_file(phased_variants)

    phased_variants.check_type('VcfFrame[Phased]')

    if imported_variants.metadata['Gene'] != phased_variants.metadata['Gene']:
        raise ValueError('Found two different genes')

    gene = imported_variants.metadata['Gene']

    if imported_variants.metadata['Assembly'] != phased_variants.metadata['Assembly']:
        raise ValueError('Found two different assemblies')

    assembly = imported_variants.metadata['Assembly']

    if imported_variants.metadata['Platform'] != phased_variants.metadata['Platform']:
        raise ValueError('Found two different platforms')

    platform = imported_variants.metadata['Platform']

    if platform in ['WGS', 'Targeted']:
        format = 'GT:AD:DP:AF'
    else:
        format = 'GT'

    vf1 = imported_variants.data.strip(format)
    vf2 = phased_variants.data.strip('GT')

    def one_row(r):
        variant = f'{r.CHROM}-{r.POS}-{r.REF}-{r.ALT}'
        s = vf1.fetch(variant)

        if s.empty:
            return r

        def one_gt(g):
            return ':'.join(g.split(':')[1:])

        s[9:] = s[9:].apply(one_gt)
        r[9:] = r[9:].str.cat(s[9:], sep=':')

        return r

    vf3 = pyvcf.VcfFrame([], vf2.df.apply(one_row, axis=1))
    vf3.df.INFO = 'Phased'
    vf3.df.FORMAT = format

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

    variant_synonyms = core.get_variant_synonyms(gene, assembly=assembly)

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

                if variant in variant_synonyms:
                    variant = variant_synonyms[variant]

                star_alleles = core.list_alleles(gene, variants=variant, assembly=assembly)

                for j in [0, 1]:
                    for star_allele in star_alleles:
                        score = 0
                        for x in anchors[sample][j]:
                            if x in core.list_variants(gene, alleles=star_allele, assembly=assembly, mode='all'):
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
    assembly='GRCh37', add_chr_prefix=False, merge=False, sv_genes=False
):
    """
    Create a BED file which contains all regions used by PyPGx.

    Parameters
    ----------
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.
    add_chr_prefix : bool, default: False
        Whether to add the 'chr' string in contig names.
    merge : bool, default: False
        Whether to merge overlapping intervals (gene names will be removed
        too).
    sv_genes : bool, default: False
        Whether to only return genes with SV.

    Returns
    -------
    fuc.api.pybed.BedFrame
        BED file.
    """
    df = core.load_gene_table()
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
    if add_chr_prefix:
        bf = bf.update_chr_prefix(mode='add')
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
        VCF file corresponding to a reference haplotype panel (zipped or
        unzipped). By default, the 1KGP panel is used.
    impute : bool, default: False
        If True, perform imputation of missing genotypes.

    Returns
    -------
    pypgx.Archive
        Archive object with the semantic type VcfFrame[Phased].
    """
    if isinstance(imported_variants, str):
        imported_variants = sdk.Archive.from_file(imported_variants)

    imported_variants.check_type('VcfFrame[Imported]')

    gene = imported_variants.metadata['Gene']
    assembly = imported_variants.metadata['Assembly']
    region = core.get_region(gene, assembly=assembly)
    beagle = f'{core.PROGRAM_PATH}/pypgx/api/beagle.28Jun21.220.jar'
    if panel is None:
        panel = f'{core.PROGRAM_PATH}/pypgx/api/1kgp/{assembly}/{gene}.vcf.gz'

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
        subprocess.run(command, check=True, stdout=subprocess.DEVNULL)
        data = pyvcf.VcfFrame.from_file(f'{t}/output.vcf.gz')
    return sdk.Archive(metadata, data)

def filter_samples(archive, samples, exclude=False):
    """
    Filter Archive for specified samples.

    Parameters
    ----------
    archive : str or pypgx.archive
        Archive file or object.
    samples : str or list
        Specify which samples should be included for analysis by providing a
        text file (.txt, .tsv, .csv, or .list) containing one sample per
        line. Alternatively, you can provide a list of samples.
    exclude : bool, default: False
        If True, exclude specified samples.

    Returns
    -------
    pypgx.Archive
        Fitlered Archive object.
    """
    if isinstance(archive, str):
        archive = sdk.Archive.from_file(archive)

    samples = common.parse_list_or_file(samples)

    if ('VcfFrame' in archive.metadata['SemanticType'] or
        'CovFrame' in archive.metadata['SemanticType']):
        data = archive.data.subset(samples, exclude=exclude)
    elif 'SampleTable' in archive.metadata['SemanticType']:
        data = archive.data.loc[samples]
    else:
        pass

    return sdk.Archive(archive.copy_metadata(), data)

def import_read_depth(
    gene, depth_of_coverage, samples=None, exclude=False
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
    samples : str or list, optional
        Specify which samples should be included for analysis by providing a
        text file (.txt, .tsv, .csv, or .list) containing one sample per
        line. Alternatively, you can provide a list of samples.
    exclude : bool, default: False
        If True, exclude specified samples.

    Returns
    -------
    pypgx.Archive
        Archive object with the semantic type CovFrame[ReadDepth].
    """
    if isinstance(depth_of_coverage, str):
        depth_of_coverage = sdk.Archive.from_file(depth_of_coverage)

    depth_of_coverage.check_type('CovFrame[DepthOfCoverage]')

    metadata = depth_of_coverage.copy_metadata()
    metadata['Gene'] = gene
    metadata['SemanticType'] = 'CovFrame[ReadDepth]'

    region = core.get_region(gene, assembly=metadata['Assembly'])

    cf = depth_of_coverage.data.update_chr_prefix(mode='remove')
    cf = cf.slice(region)

    if samples is not None:
        samples = common.parse_list_or_file(samples)
        cf = cf.subset(samples, exclude=exclude)

    return sdk.Archive(metadata, cf)

def import_variants(
    gene, vcf, assembly='GRCh37', platform='WGS', samples=None, exclude=False
):
    """
    Import variant (SNV/indel) data for the target gene.

    The method will first slice input VCF for the target gene and then assess
    whether every genotype call in the sliced VCF is haplotype phased. It
    will return an archive object with the semantic type
    VcfFrame[Consolidated] if the VCF is fully phased or otherwise
    VcfFrame[Imported].

    Parameters
    ----------
    gene : str
        Target gene.
    vcf : str or fuc.api.pyvcf.VcfFrame
        Input VCF file must be already BGZF compressed (.gz) and indexed
        (.tbi) to allow random access. Alternatively, you can provide a
        VcfFrame object.
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.
    platform : {'WGS', 'Targeted', 'Chrip'}, default: 'WGS'
        Genotyping platform.
    samples : str or list, optional
        Specify which samples should be included for analysis by providing a
        text file (.txt, .tsv, .csv, or .list) containing one sample per
        line. Alternatively, you can provide a list of samples.
    exclude : bool, default: False
        If True, exclude specified samples.

    Returns
    -------
    pypgx.Archive
        Archive object with the semantic type VcfFrame[Imported] or
        VcfFrame[Consolidated].
    """
    region = core.get_region(gene, assembly=assembly)

    if isinstance(vcf, str):
        vf = pyvcf.VcfFrame.from_file(vcf, regions=region)
    else:
        vf = vcf.slice(region)

    vf = vf.update_chr_prefix(mode='remove')
    vf = vf.strip('GT:AD:DP')
    vf = vf.add_af()

    if samples is not None:
        samples = common.parse_list_or_file(samples)
        vf = vf.subset(samples, exclude=exclude)

    if vf.phased:
        semantic_type = 'VcfFrame[Consolidated]'
    else:
        vf = vf.unphase()
        semantic_type = 'VcfFrame[Imported]'

    metadata = {
        'Platform': platform,
        'Gene': gene,
        'Assembly': assembly,
        'SemanticType': semantic_type,
    }

    return sdk.Archive(metadata, vf)

def predict_alleles(consolidated_variants):
    """
    Predict candidate star alleles based on observed SNVs and indels.

    Parameters
    ----------
    consolidated_variants : str or pypgx.Archive
        Archive file or object with the semantic type VcfFrame[Consolidated].

    Returns
    -------
    pypgx.Archive
        Archive object with the semantic type SampleTable[Alleles].
    """
    if isinstance(consolidated_variants, str):
        consolidated_variants = sdk.Archive.from_file(consolidated_variants)

    consolidated_variants.check_type('VcfFrame[Consolidated]')

    gene = consolidated_variants.metadata['Gene']
    assembly = consolidated_variants.metadata['Assembly']

    definition_table = core.build_definition_table(gene, assembly)
    ref_allele = core.get_ref_allele(gene, assembly)
    default_allele = core.get_default_allele(gene, assembly)
    defining_variants = core.list_variants(gene, assembly=assembly)
    variant_synonyms = core.get_variant_synonyms(gene, assembly=assembly)

    reformatted_variants = {}

    for x in consolidated_variants.data.variants():
        if x in variant_synonyms:
            y = variant_synonyms[x]
            if y in reformatted_variants:
                raise ValueError('Multiple variant synonyms detected')
            reformatted_variants[y] = x

    star_alleles = {}

    for allele in definition_table.samples:
        df = definition_table.df[definition_table.df[allele] == '1']
        star_alleles[allele] = set(df.apply(lambda r: f'{r.CHROM}-{r.POS}-{r.REF}-{r.ALT}', axis=1))

    samples = {}

    def one_haplotype(observed):
        """
        Call candidate alleles for haplotype.
        """
        candidates = []
        for allele, variants in star_alleles.items():
            if variants.issubset(observed):
                candidates.append(allele)
        candidates = core.collapse_alleles(gene, candidates, assembly=assembly)
        if ref_allele != default_allele and ref_allele not in candidates and default_allele not in candidates:
            candidates.append(default_allele)
        if not candidates:
            candidates.append(default_allele)
        candidates = core.sort_alleles(candidates, by='priority', gene=gene, assembly=assembly)
        return candidates

    def one_row(r, sample, i):
        gt = r[sample].split(':')[0]
        if '.' in gt:
            return ''
        j = int(gt.split('|')[i])
        if j == 0:
            return ''
        alt = r.ALT.split(',')[j-1]
        variant = f'{r.CHROM}-{r.POS}-{r.REF}-{alt}'
        if variant in variant_synonyms:
            variant = variant_synonyms[variant]
        if variant not in defining_variants:
            return ''
        return variant

    for sample in consolidated_variants.data.samples:
        results = []
        alt_phase = []
        all_alleles = []

        for i in [0, 1, 2]:
            if i == 2:
                candidates = one_haplotype(set(alt_phase))
                candidates = [x for x in candidates if x not in all_alleles]
                all_alleles += [x for x in candidates if x not in all_alleles]
                all_alleles = core.sort_alleles(all_alleles, by='priority', gene=gene, assembly=assembly)
            else:
                observed = consolidated_variants.data.df.apply(one_row, args=(sample, i), axis=1)
                observed = [x for x in observed if x]
                alt_phase += [x for x in observed if x not in alt_phase]
                candidates = one_haplotype(observed)
                all_alleles += [x for x in candidates if x not in all_alleles]

            results.append(';'.join(candidates) + ';')

        af_list = []

        for allele in all_alleles:
            if allele == default_allele:
                af_list.append(f'{allele}:default')
            else:
                variants = ','.join(star_alleles[allele])
                fractions = ','.join([str(consolidated_variants.data.get_af(sample, reformatted_variants[x])) if x in reformatted_variants else str(consolidated_variants.data.get_af(sample, x)) for x in star_alleles[allele]])
                af_list.append(f'{allele}:{variants}:{fractions}')

        results.append(';'.join(af_list) + ';')
        samples[sample] = results

    data = pd.DataFrame(samples).T
    data.columns = ['Haplotype1', 'Haplotype2', 'AlternativePhase', 'VariantData']
    metadata = consolidated_variants.copy_metadata()
    metadata['SemanticType'] = 'SampleTable[Alleles]'

    return sdk.Archive(metadata, data)

def predict_cnv(copy_number, cnv_caller=None):
    """
    Predict CNV for the target gene based on copy number data.

    Genomic positions that are missing copy number, because for example the
    input data is targeted sequencing, will be imputed with forward filling.

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

    copy_number.check_type('CovFrame[CopyNumber]')

    gene = copy_number.metadata['Gene']
    assembly = copy_number.metadata['Assembly']
    model_file = f'{core.PROGRAM_PATH}/pypgx/api/cnv/{assembly}/{gene}.zip'

    if cnv_caller is None:
        cnv_caller = sdk.Archive.from_file(model_file)
    else:
        if isinstance(cnv_caller, str):
            cnv_caller = sdk.Archive.from_file(cnv_caller)

        cnv_caller.check_type('Model[CNV]')

    copy_number = _process_copy_number(copy_number)

    df = copy_number.data.df.iloc[:, 2:]
    X = df.T.to_numpy()
    predictions = cnv_caller.data.predict(X)
    df = core.load_cnv_table()
    df = df[df.Gene == copy_number.metadata['Gene']]
    cnvs = dict(zip(df.Code, df.Name))
    predictions = [cnvs[x] for x in predictions]
    metadata = copy_number.copy_metadata()
    metadata['SemanticType'] = 'SampleTable[CNVCalls]'
    data = pd.DataFrame({'CNV': predictions})
    data.index = copy_number.data.samples
    return sdk.Archive(metadata, data)

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
        merge=True, sv_genes=True, assembly=assembly,
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
            bf = bf.update_chr_prefix(mode='add')
        else:
            bf = bf.update_chr_prefix(mode='remove')
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

    cnv_caller.check_type('Model[CNV]')

    if isinstance(copy_number, str):
        copy_number = sdk.Archive.from_file(copy_number)

    copy_number.check_type('CovFrame[CopyNumber]')

    if isinstance(cnv_calls, str):
        cnv_calls = sdk.Archive.from_file(cnv_calls)

    cnv_calls.check_type('SampleTable[CNVCalls]')

    sdk.compare_metadata('Gene', cnv_caller, copy_number, cnv_calls)
    sdk.compare_metadata('Assembly', cnv_caller, copy_number, cnv_calls)

    copy_number = _process_copy_number(copy_number)

    cnv_table = core.load_cnv_table()
    cnv_table = cnv_table[cnv_table.Gene == copy_number.metadata['Gene']]
    name2code = dict(zip(cnv_table.Name, cnv_table.Code))
    code2name = dict(zip(cnv_table.Code, cnv_table.Name))

    cnv_calls.data['Code'] = cnv_calls.data.apply(lambda r: name2code[r.CNV], axis=1)
    columns = ['Chromosome', 'Position'] + cnv_calls.data.Sample.to_list()
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
    calls using the one-vs-rest strategy.

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

    copy_number.check_type('CovFrame[CopyNumber]')

    if isinstance(cnv_calls, str):
        cnv_calls = sdk.Archive.from_file(cnv_calls)

    cnv_calls.check_type('SampleTable[CNVCalls]')

    sdk.compare_metadata('Gene', copy_number, cnv_calls)
    sdk.compare_metadata('Assembly', copy_number, cnv_calls)

    copy_number = _process_copy_number(copy_number)

    cnv_table = core.load_cnv_table()
    cnv_table = cnv_table[cnv_table.Gene == copy_number.metadata['Gene']]
    name2code = dict(zip(cnv_table.Name, cnv_table.Code))
    code2name = dict(zip(cnv_table.Code, cnv_table.Name))
    cnv_calls.data['Code'] = cnv_calls.data.apply(lambda r: name2code[r.CNV], axis=1)
    columns = ['Chromosome', 'Position'] + cnv_calls.data.Sample.to_list()
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
