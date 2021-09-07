from . import utils
from .. import sdk

from fuc import pyvcf, pycov, common
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

###################
# Private methods #
###################

def _plot_exons(gene, assembly, ax):
    region = utils.get_region(gene, assembly=assembly)
    chrom, start, end = common.parse_region(region)
    df = utils.load_gene_table()
    starts1 = [int(x) for x in df[df.Gene == gene][f'{assembly}ExonStarts'].values[0].strip(',').split(',')]
    ends1 = [int(x) for x in df[df.Gene == gene][f'{assembly}ExonStarts'].values[0].strip(',').split(',')]
    paralog = utils.get_paralog(gene)
    if paralog:
        starts2 = [int(x) for x in df[df.Gene == paralog][f'{assembly}ExonStarts'].values[0].strip(',').split(',')]
        ends2 = [int(x) for x in df[df.Gene == paralog][f'{assembly}ExonStarts'].values[0].strip(',').split(',')]
    common.plot_exons(
        starts1, ends1, ax=ax, name=gene, fontsize=20
    )
    if paralog:
        common.plot_exons(
            starts2, ends2, ax=ax, name=paralog, fontsize=20
    )
    ax.set_xlim([start, end])
    ax.axis('off')

##################
# Public methods #
##################

def plot_bam_copy_number(
    result, path=None, samples=None, ymin=None, ymax=None
):
    """
    Plot copy number profile with BAM data.

    Parameters
    ----------
    result : pypgx.Archive or str
        Archive file with the semantic type CovFrame[CopyNumber].
    path : str, optional
        Create plots in this directory.
    samples : list, optional
        Create plots only for these samples.
    ymin : float, optional
        Y-axis bottom.
    ymax : float, optional
        Y-axis top.
    """
    if isinstance(result, str):
        result = sdk.Archive.from_file(result)

    if result.metadata['SemanticType'] != 'CovFrame[CopyNumber]':
        raise ValueError('Incorrect semantic type')

    if samples is None:
        samples = result.data.samples

    with sns.axes_style('darkgrid'):
        for sample in samples:

            fig, [ax1, ax2] = plt.subplots(2, 1, figsize=(18, 12), gridspec_kw={'height_ratios': [1, 10]})

            _plot_exons(result.metadata['Gene'], result.metadata['Assembly'], ax1)
            result.data.plot_region(sample, ax=ax2, legend=False)

            ax2.set_ylim([ymin, ymax])
            ax2.set_xlabel('Coordinate', fontsize=25)
            ax2.set_ylabel('Copy number', fontsize=25)
            ax2.tick_params(axis='both', which='major', labelsize=20)

            if path is None:
                output = f'{sample}.png'
            else:
                output = f'{path}/{sample}.png'

            plt.tight_layout()
            fig.savefig(output)
            plt.close()

def plot_bam_read_depth(
    result, path=None, samples=None, ymin=None, ymax=None
):
    """
    Plot copy number profile with BAM data.

    Parameters
    ----------
    result : pypgx.Archive or str
        Archive file with the semantic type CovFrame[ReadDepth].
    path : str, optional
        Create plots in this directory.
    samples : list, optional
        Create plots only for these samples.
    ymin : float, optional
        Y-axis bottom.
    ymax : float, optional
        Y-axis top.
    """

    if isinstance(result, str):
        result = sdk.Archive.from_file(result)

    if result.metadata['SemanticType'] != 'CovFrame[ReadDepth]':
        raise ValueError('Incorrect semantic type')

    if samples is None:
        samples = result.data.samples

    with sns.axes_style('darkgrid'):
        for sample in samples:

            fig, [ax1, ax2] = plt.subplots(2, 1, figsize=(18, 12), gridspec_kw={'height_ratios': [1, 10]})

            _plot_exons(result.metadata['Gene'], result.metadata['Assembly'], ax1)

            result.data.plot_region(sample, ax=ax2, legend=False)

            ax2.set_ylim([ymin, ymax])
            ax2.set_xlabel('Coordinate', fontsize=25)
            ax2.set_ylabel('Read depth', fontsize=25)
            ax2.tick_params(axis='both', which='major', labelsize=20)

            if path is None:
                output = f'{sample}.png'
            else:
                output = f'{path}/{sample}.png'

            plt.tight_layout()
            fig.savefig(output)
            plt.close()

def plot_vcf_allele_fraction(
    gene, vcf, assembly='GRCh37', path=None, samples=None, ymin=None,
    ymax=None
):
    """
    Plot allele fraction profile with VCF data.

    Parameters
    ----------
    gene : str
        Target gene.
    vcf : str
        VCF file.
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.
    path : str, optional
        Create plots in this directory.
    samples : list, optional
        Create plots only for these samples.
    ymin : float, optional
        Y-axis bottom.
    ymax : float, optional
        Y-axis top.
    """
    vf = pyvcf.VcfFrame.from_file(vcf)

    region = utils.get_region(gene, assembly=assembly)
    chrom, start, end = common.parse_region(region)

    if samples is None:
        samples = cf.samples

    with sns.axes_style('darkgrid'):
        for sample in samples:

            fig, [ax1, ax2] = plt.subplots(2, 1, figsize=(18, 12), gridspec_kw={'height_ratios': [1, 10]})

            _plot_exons(gene, assembly, ax1)

            vf.plot_region(sample, region=region, ax=ax2, k='#AD_FRAC_REF', label='REF')
            vf.plot_region(sample, region=region, ax=ax2, k='#AD_FRAC_ALT', label='ALT')

            if path is None:
                output = f'{sample}.png'
            else:
                output = f'{path}/{sample}.png'

            ax2.set_ylim([ymin, ymax])
            ax2.set_xlabel(f'Chromosome {chrom}', fontsize=25)
            ax2.set_ylabel('Allele fraction', fontsize=25)
            ax2.tick_params(axis='both', which='major', labelsize=20)

            plt.tight_layout()
            fig.savefig(output)
            plt.close()

def plot_vcf_read_depth(
    gene, vcf, assembly='GRCh37', path=None, samples=None, ymin=None,
    ymax=None
):
    """
    Plot read depth profile with VCF data.

    Parameters
    ----------
    gene : str
        Target gene.
    vcf : str
        VCF file.
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.
    path : str, optional
        Create plots in this directory.
    samples : list, optional
        Create plots only for these samples.
    ymin : float, optional
        Y-axis bottom.
    ymax : float, optional
        Y-axis top.
    """
    vf = pyvcf.VcfFrame.from_file(vcf)

    region = utils.get_region(gene, assembly=assembly)
    chrom, start, end = common.parse_region(region)

    if samples is None:
        samples = cf.samples

    with sns.axes_style('darkgrid'):
        for sample in samples:

            fig, [ax1, ax2] = plt.subplots(2, 1, figsize=(18, 12), gridspec_kw={'height_ratios': [1, 10]})

            _plot_exons(gene, assembly, ax1)

            vf.plot_region(sample, region=region, ax=ax2, label='ALT')

            if path is None:
                output = f'{sample}.png'
            else:
                output = f'{path}/{sample}.png'

            ax2.set_ylim([ymin, ymax])
            ax2.set_xlabel(f'Chromosome {chrom}', fontsize=25)
            ax2.set_ylabel('Read depth', fontsize=25)
            ax2.tick_params(axis='both', which='major', labelsize=20)

            plt.tight_layout()
            fig.savefig(output)
            plt.close()
