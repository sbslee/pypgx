"""
The utils submodule contains main actions of PyPGx.
"""

import pkgutil
from io import BytesIO
import tempfile
import zipfile
import subprocess
import os
import sys
import pickle
import warnings

from . import core
from .. import sdk

import numpy as np
import pandas as pd
import pysam
from fuc import pybam, pyvcf, pycov, common, pybed
from sklearn import metrics
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import SVC
from scipy.ndimage import median_filter

###################
# Private methods #
###################

def _phase_extension(vf, gene, assembly):
    """
    Apply the phase-extension algorithm.

    Anchor variants are those variants that have been haplotype phased
    using a reliable method (e.g. statistical haplotype phasing and
    read-backed phasing) and are later used by the phase-extension
    algorithm (PE). Basically, PE determines the most likely haplotype
    phase of the remaining unphased variants using anchor variants.
    For each unphased variant, PE first finds all star alleles carrying the
    variant and then counts how many anchor variants per haplotype are
    overlapped to each of the star alleles. For example, if the second
    haplotype's anchor variants (i.e. variants with '0|1') were found to
    have the most overlapping with the *2 allele, then PE will assign the
    phase of the variant of interest to '0|1'.
    """
    anchors = {}

    for i, r in vf.df.iterrows():
        for allele in r.ALT.split(','):
            variant = f'{r.CHROM}-{r.POS}-{r.REF}-{allele}'
            for sample in vf.samples:
                if sample not in anchors:
                    anchors[sample] = [[], []]
                gt = r[sample].split(':')[0]
                if '|' not in gt:
                    continue
                gt = gt.split('|')
                if gt[0] != '0':
                    anchors[sample][0].append(variant)
                if gt[1] != '0':
                    anchors[sample][1].append(variant)

    variant_synonyms = core.get_variant_synonyms(gene, assembly=assembly)

    def one_row(r):
        if pyvcf.row_phased(r):
            return r

        r.FORMAT += ':PE'

        for sample in vf.samples:
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

    return pyvcf.VcfFrame([], vf.df.apply(one_row, axis=1))

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
    Call phenotypes for target gene.

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
    Calculate concordance between two genotype results.

    Only samples that appear in both genotype results will be used to
    calculate concordance for genotype calls as well as CNV calls.

    Parameters
    ----------
    first : str or pypgx.Archive
        First archive file or object with the semantic type
        SampleTable[Results].
    second : str or pypgx.Archive
        Second archive file or object with the semantic type
        SampleTable[Results].
    verbose : bool, default: False
        If True, print the verbose version of output, including discordant
        calls.

    Examples
    --------

    >>> import pypgx
    >>> pypgx.compare_genotypes('results-1.zip', 'results-2.zip')
    # Genotype
    Total: 100
    Compared: 100
    Concordance: 1.000 (100/100)
    # CNV
    Total: 100
    Compared: 100
    Concordance: 1.000 (100/100)
    """
    if isinstance(first, str):
        first = sdk.Archive.from_file(first)

    first.check_type('SampleTable[Results]')

    if isinstance(second, str):
        second = sdk.Archive.from_file(second)

    second.check_type('SampleTable[Results]')

    def show_comparison(col):
        df = pd.concat([first.data[col], second.data[col]], axis=1)
        print(f'# {col}')
        print(f'Total: {df.shape[0]}')
        df.columns = ['First', 'Second']
        df = df.dropna()
        print(f'Compared: {df.shape[0]}')
        df['Concordant'] = df.First == df.Second
        if len(df.Concordant):
            print(f'Concordance: {sum(df.Concordant)/len(df.Concordant):.3f} ({sum(df.Concordant)}/{len(df.Concordant)})')
        else:
            print('Concordance: N/A')
        if verbose:
            print('Discordant:')
            if df.Concordant.all():
                print('None')
            else:
                print(df[~df.Concordant])

    for col in ['Genotype', 'CNV']:
        show_comparison(col)

def compute_control_statistics(
    gene, bams, assembly='GRCh37', bed=None
):
    """
    Compute summary statistics for control gene from BAM files.

    Note that for the arguments ``gene`` and ``bed``, the 'chr' prefix in
    contig names (e.g. 'chr1' vs. '1') will be automatically added or removed
    as necessary to match the input BAM's contig names.

    Parameters
    ----------
    gene : str
        Control gene (recommended choices: 'EGFR', 'RYR1', 'VDR').
        Alternatively, you can provide a custom region (format:
        chrom:start-end).
    bams : str or list
        One or more input BAM files. Alternatively, you can provide a text
        file (.txt, .tsv, .csv, or .list) containing one BAM file per line.
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.
    bed : str, optional
        By default, the input data is assumed to be WGS. If it's targeted
        sequencing, you must provide a BED file to indicate probed regions.

    Returns
    -------
    pypgx.Archive
        Archive object with the semantic type SampleTable[Statistics].
    """
    gene_table = core.load_gene_table()

    if gene in core.list_genes(mode='all'):
        region = gene_table[gene_table.Gene == gene][f'{assembly}Region'].values[0]
    else:
        region = gene

    cf = pycov.CovFrame.from_bam(bams, regions=region, zero=False)

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
    control_statistics : str or pypgx.Archive
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
    gene, bams, assembly='GRCh37', bed=None
):
    """
    Compute read depth for target gene from BAM files.

    By default, the input data is assumed to be WGS. If it's targeted
    sequencing, you must provide a BED file with ``bed`` to indicate
    probed regions.

    Parameters
    ----------
    gene : str
        Target gene.
    bams : str or list
        One or more input BAM files. Alternatively, you can provide a text
        file (.txt, .tsv, .csv, or .list) containing one BAM file per line.
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

    region = core.get_region(gene, assembly=assembly)

    data = pycov.CovFrame.from_bam(bams, regions=region, zero=True)

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

    # For every variant in VcfFrame[Phased] (e.g. '0|1'), find and append its
    # accompanying data from VcfFrame[Imported] (e.g. '0|1:15,15:30:0.5').
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

    # Remove variants that are in both VcfFrame[Imported] and
    # VcfFrame[Phased]. Append remaining unphased variants to
    # VcfFrame[Phased].
    vf4 = vf1.filter_vcf(vf2, opposite=True)
    vf5 = pyvcf.VcfFrame([], pd.concat([vf3.df, vf4.df])).sort()

    vf6 = _phase_extension(vf5, gene, assembly)

    metadata = phased_variants.copy_metadata()
    metadata['SemanticType'] = 'VcfFrame[Consolidated]'

    return sdk.Archive(metadata, vf6)

def create_input_vcf(
    vcf, fasta, bams, assembly='GRCh37', genes=None, exclude=False,
    dir_path=None, max_depth=250
):
    """
    Call SNVs/indels from BAM files for all target genes.

    To save computing resources, this method will call variants only for
    target genes whose at least one star allele is defined by SNVs/indels.
    Therefore, variants will not be called for target genes that have star
    alleles defined only by structural variation (e.g. UGT2B17).

    Parameters
    ----------
    vcf : str
        Output VCF file. It must have .vcf.gz as suffix.
    fasta : str
        Reference FASTA file.
    bams : str or list
        One or more input BAM files. Alternatively, you can provide a text
        file (.txt, .tsv, .csv, or .list) containing one BAM file per line.
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.
    genes : list, optional
        List of genes to include.
    exclude : bool, default: False
        Exclude specified genes. Ignored when ``genes=None``.
    dir_path : str, optional
        By default, intermediate files (likelihoods.bcf, calls.bcf, and
        calls.normalized.bcf) will be stored in a temporary directory, which
        is automatically deleted after creating final VCF. If you provide a
        directory path, intermediate files will be stored there.
    max_depth : int, default: 250
        At a position, read maximally this number of reads per input file.
        If your input data is from WGS (e.g. 30X), you don't need to change
        this option. However, if it's from targeted sequencing with
        ultra-deep coverage (e.g. 500X), then you need to increase the
        maximum depth.
    """
    if not vcf.endswith('.vcf.gz'):
        raise ValueError(f"VCF file must have .vcf.gz as suffix: {vcf}")
    vcf = vcf.replace('.vcf.gz', '.vcf')
    bf = create_regions_bed(merge=True, assembly=assembly, var_genes=True,
        genes=genes, exclude=exclude)
    pyvcf.call(fasta=fasta, bams=bams, regions=bf, path=vcf, gap_frac=0,
        dir_path=dir_path, group_samples='-', max_depth=max_depth)
    pysam.tabix_index(vcf, preset='vcf', force=True)

def create_regions_bed(
    assembly='GRCh37', add_chr_prefix=False, merge=False, target_genes=False,
    sv_genes=False, var_genes=False, genes=None, exclude=False
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
    target_genes : bool, default: False
        Whether to only return target genes, excluding control genes and
        paralogs.
    sv_genes : bool, default: False
        Whether to only return target genes whose at least one star allele is
        defined by structural variation.
    var_genes : bool, default: False
        Whether to only return target genes whose at least one star allele is
        defined by SNVs/indels.
    genes : list, optional
        List of genes to include.
    exclude : bool, default: False
        Exclude specified genes. Ignored when ``genes=None``.

    Returns
    -------
    fuc.api.pybed.BedFrame
        BED file.

    Examples
    --------

    >>> import pypgx
    >>> bf = pypgx.create_regions_bed()
    >>> bf.gr.df.head()
      Chromosome      Start        End     Name
    0          1  201005639  201084694  CACNA1S
    1          1   60355979   60395470   CYP2J2
    2          1   47391859   47410148  CYP4A11
    3          1   47600112   47618399  CYP4A22
    4          1   47261669   47288021   CYP4B1
    >>> bf = pypgx.create_regions_bed(assembly='GRCh38')
    >>> bf.gr.df.head()
      Chromosome      Start        End     Name
    0          1  201036511  201115426  CACNA1S
    1          1   59890307   59929773   CYP2J2
    2          1   46926187   46944476  CYP4A11
    3          1   47134440   47152727  CYP4A22
    4          1   46796045   46822413   CYP4B1
    >>> bf = pypgx.create_regions_bed(add_chr_prefix=True)
    >>> bf.gr.df.head()
      Chromosome      Start        End     Name
    0       chr1  201005639  201084694  CACNA1S
    1       chr1   60355979   60395470   CYP2J2
    2       chr1   47391859   47410148  CYP4A11
    3       chr1   47600112   47618399  CYP4A22
    4       chr1   47261669   47288021   CYP4B1
    >>> bf = pypgx.create_regions_bed(merge=True)
    >>> bf.gr.df.head()
      Chromosome     Start       End
    0          1  47261669  47288021
    1          1  47391859  47410148
    2          1  47600112  47618399
    3          1  60355979  60395470
    4          1  97540298  98389615

    """
    df = core.load_gene_table()
    if genes is not None:
        if exclude:
            df = df[~df.Gene.isin(genes)]
        else:
            df = df[df.Gene.isin(genes)]
    if target_genes:
        df = df[df.Target]
    if sv_genes:
        df = df[df.SV]
    if var_genes:
        df = df[df.Variants]
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
        Archive file or object with the semantic type VcfFrame[Imported]. The
        'chr' prefix in contig names (e.g. 'chr1' vs. '1') will be
        automatically added or removed as necessary to match the reference
        VCF's contig names.
    panel : str, optional
        VCF file corresponding to a reference haplotype panel (compressed or
        uncompressed). By default, the 1KGP panel in the ``~/pypgx-bundle``
        directory will be used.
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

    metadata = imported_variants.copy_metadata()
    metadata['SemanticType'] = 'VcfFrame[Phased]'
    metadata['Program'] = 'Beagle'

    if panel is None:
        home = os.path.expanduser('~')
        panel = f'{home}/pypgx-bundle/1kgp/{assembly}/{gene}.vcf.gz'

    has_chr_prefix = pyvcf.has_chr_prefix(panel)

    if has_chr_prefix:
        vf1 = imported_variants.data.update_chr_prefix('add')
        region = common.update_chr_prefix(region, 'add')
    else:
        vf1 = imported_variants.data.copy()

    if vf1.empty:
        return sdk.Archive(metadata, vf1)

    with tempfile.TemporaryDirectory() as t:
        vf1.to_file(f'{t}/input.vcf')
        command = [
            'java', '-Xmx2g', '-jar', beagle,
            f'gt={t}/input.vcf',
            f'chrom={region}',
            f'ref={panel}',
            f'out={t}/output',
            f'impute={str(impute).lower()}'
        ]
        subprocess.run(command, check=True, stdout=subprocess.DEVNULL)
        vf2 = pyvcf.VcfFrame.from_file(f'{t}/output.vcf.gz')

    if has_chr_prefix:
        vf2 = vf2.update_chr_prefix('remove')

    return sdk.Archive(metadata, vf2)

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
        if exclude:
            samples = [x for x in archive.data.index.to_list()
                if x not in samples]
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
    Import SNV/indel data for target gene.

    The method will slice the input VCF for the target gene to create an
    archive object with the semantic type VcfFrame[Imported] or
    VcfFrame[Consolidated].

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
    platform : {'WGS', 'Targeted', 'Chip', 'LongRead'}, default: 'WGS'
        Genotyping platform used. When the platform is 'WGS', 'Targeted', or
        'Chip', the method will assess whether every genotype call in the
        sliced VCF is haplotype phased (e.g. '0|1'). If the sliced VCF is
        fully phased, the method will return VcfFrame[Consolidated] or
        otherwise VcfFrame[Imported]. When the platform is 'LongRead', the
        method will return VcfFrame[Consolidated] after applying the
        phase-extension algorithm to estimate haplotype phase of any variants
        that could not be resolved by read-backed phasing.
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

    if platform == 'LongRead':
        vf = _phase_extension(vf, gene, assembly)
        semantic_type = 'VcfFrame[Consolidated]'
    else:
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

    for x in consolidated_variants.data.to_variants():
        if x in variant_synonyms:
            y = variant_synonyms[x]
            if y in reformatted_variants:
                warnings.warn(f"Multiple variant synonyms detected for {y}: PyPGx will report information for {x}")
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
    Predict CNV from copy number data for target gene.

    Genomic positions that are missing copy number because, for example, the
    input data is targeted sequencing will be imputed with forward filling.

    Parameters
    ----------
    copy_number : str or pypgx.Archive
        Archive file or object with the semantic type CovFrame[CopyNumber].
    cnv_caller : str or pypgx.Archive, optional
        Archive file or object with the semantic type Model[CNV]. By default,
        a pre-trained CNV caller in the ``~/pypgx-bundle`` directory will be
        used.

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
    home = os.path.expanduser('~')
    model_file = f'{home}/pypgx-bundle/cnv/{assembly}/{gene}.zip'

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
    cnv_table = core.load_cnv_table()
    cnv_table = cnv_table[cnv_table.Gene == copy_number.metadata['Gene']]
    code2name = dict(zip(list(range(len(cnv_table.Name))), cnv_table.Name))
    predictions = [code2name[x] for x in predictions]
    metadata = copy_number.copy_metadata()
    metadata['SemanticType'] = 'SampleTable[CNVCalls]'
    data = pd.DataFrame({'CNV': predictions})
    data.index = copy_number.data.samples

    return sdk.Archive(metadata, data)

def prepare_depth_of_coverage(
    bams, assembly='GRCh37', bed=None, genes=None, exclude=False
):
    """
    Prepare a depth of coverage file for all target genes with SV from BAM
    files.

    To save computing resources, this method will count read depth only for
    target genes whose at least one star allele is defined by structural
    variation. Therefore, read depth will not be computed for target genes
    that have star alleles defined only by SNVs/indels (e.g. CYP3A5).

    Parameters
    ----------
    bams : str or list
        One or more input BAM files. Alternatively, you can provide a text
        file (.txt, .tsv, .csv, or .list) containing one BAM file per line.
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.
    bed : str, optional
        By default, the input data is assumed to be WGS. If it's targeted
        sequencing, you must provide a BED file to indicate probed regions.
        Note that the 'chr' prefix in contig names (e.g. 'chr1' vs. '1') will
        be automatically added or removed as necessary to match the input
        BAM's contig names.
    genes : list, optional
        List of genes to include.
    exclude : bool, default: False
        Exclude specified genes. Ignored when ``genes=None``.

    Returns
    -------
    pypgx.Archive
        Archive object with the semantic type CovFrame[DepthOfCoverage].
    """
    metadata = {
        'Assembly': assembly,
        'SemanticType': 'CovFrame[DepthOfCoverage]',
    }

    regions = create_regions_bed(
        merge=True, sv_genes=True, assembly=assembly, genes=genes,
        exclude=exclude
    ).to_regions()

    cf = pycov.CovFrame.from_bam(bams, regions=regions, zero=True)

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

def print_data(input):
    """
    Print the main data of specified archive.

    Parameters
    ----------
    input : pypgx.Archive
        Archive file.
    """
    archive = sdk.Archive.from_file(input)
    if 'SampleTable' in archive.type:
        data = archive.data.to_csv(sep='\t')
    elif 'CovFrame' in archive.type:
        data = archive.data.to_string()
    elif 'VcfFrame' in archive.type:
        data = archive.data.to_string()
    else:
        raise ValueError(f"Data cannot be printed for {archive.type}")

    # https://docs.python.org/3/library/signal.html#note-on-sigpipe
    try:
        print(data, end='')
    except BrokenPipeError:
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        sys.exit(1)

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

def slice_bam(
    input, output, assembly='GRCh37', genes=None, exclude=False
):
    """
    Slice BAM file for all genes used by PyPGx.

    Parameters
    ----------
    input
        Input BAM file. It must be already indexed to allow random access.
    output : str
        Output BAM file.
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.
    genes : list, optional
        List of genes to include.
    exclude : bool, default: False
        Exclude specified genes. Ignored when ``genes=None``.
    """
    bf = create_regions_bed(merge=True, assembly=assembly,
        genes=genes, exclude=exclude)
    pybam.slice(input, bf, path=output)

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
        Write the confusion matrix as a CSV file where rows indicate actual
        class and columns indicate prediction class.
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
    code = list(range(len(cnv_table.Name)))
    code2name = dict(zip(code, cnv_table.Name))
    name2code = dict(zip(cnv_table.Name, code))
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
        Write the confusion matrix as a CSV file where rows indicate actual
        class and columns indicate prediction class.

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
    code = list(range(len(cnv_table.Name)))
    code2name = dict(zip(code, cnv_table.Name))
    name2code = dict(zip(cnv_table.Name, code))
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
