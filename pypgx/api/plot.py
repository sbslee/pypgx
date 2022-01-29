"""
The plot submodule is used to plot various kinds of profiles such as read
depth, copy number, and allele fraction.
"""

from . import utils, core
from .. import sdk

from fuc import pyvcf, pycov, common
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

###################
# Private methods #
###################

def _plot_exons(gene, assembly, ax, fontsize=25):
    """Plot a gene model."""
    region = core.get_region(gene, assembly=assembly)
    chrom, start, end = common.parse_region(region)
    df = core.load_gene_table()
    starts1 = core.get_exon_starts(gene, assembly=assembly)
    ends1 = core.get_exon_ends(gene, assembly=assembly)
    paralog = core.get_paralog(gene)
    strand = core.get_strand(gene)
    name1 = f'{gene} ({strand})'
    if paralog:
        starts2 = core.get_exon_starts(paralog, assembly=assembly)
        ends2 = core.get_exon_ends(paralog, assembly=assembly)
        name2 = f'{paralog} ({strand})'
    common.plot_exons(
        starts1, ends1, ax=ax, name=name1, fontsize=fontsize, offset=2
    )
    if paralog:
        common.plot_exons(
            starts2, ends2, ax=ax, name=name2, fontsize=fontsize, offset=2
    )
    ax.set_ylim([-1.5, 1.5])
    ax.set_xlim([start, end])
    ax.axis('off')

def _plot_bam_copy_number_one(
    ax1, ax2, sample, copy_number, gene, assembly, processed_copy_number,
    ymin, ymax, fontsize
):
    region = core.get_region(gene, assembly=assembly)
    chrom, start, end = common.parse_region(region)

    _plot_exons(gene, assembly, ax1, fontsize=fontsize)

    copy_number.data.plot_region(sample, ax=ax2, legend=False)

    if processed_copy_number is not None:
        processed_copy_number.data.plot_region(sample,
            ax=ax2, legend=False)

    ax2.set_xlim([start, end])
    ax2.locator_params(axis='x', nbins=4)
    ax2.set_ylim([ymin, ymax])
    ax2.set_xlabel(f'Coordinate in chr{chrom} (Mb)', fontsize=fontsize)
    ax2.set_ylabel('Copy number', fontsize=fontsize)
    ax2.tick_params(axis='both', which='major', labelsize=fontsize)
    ax2.ticklabel_format(axis='x', useOffset=False, scilimits=(6, 6))
    ax2.xaxis.get_offset_text().set_fontsize(fontsize)

    return ax1, ax2

def _plot_vcf_allele_fraction_one(
    ax1, ax2, sample, imported_variants, gene, assembly, fontsize
):
    region = core.get_region(gene, assembly=assembly)
    chrom, start, end = common.parse_region(region)

    _plot_exons(gene, assembly, ax1, fontsize=fontsize)

    imported_variants.data.plot_region(sample, ax=ax2, k='#AD_FRAC_REF', label='REF')
    imported_variants.data.plot_region(sample, ax=ax2, k='#AD_FRAC_ALT', label='ALT')

    ax2.set_xlim([start, end])
    ax2.locator_params(axis='x', nbins=4)
    ax2.set_ylim([-0.05, 1.05])
    ax2.set_xlabel(f'Coordinate in chr{chrom} (Mb)', fontsize=fontsize)
    ax2.set_ylabel('Allele fraction', fontsize=fontsize)
    ax2.tick_params(axis='both', which='major', labelsize=fontsize)
    ax2.ticklabel_format(axis='x', useOffset=False, scilimits=(6, 6))
    ax2.xaxis.get_offset_text().set_fontsize(fontsize)

    return ax1, ax2

##################
# Public methods #
##################

def plot_bam_copy_number(
    copy_number, fitted=False, path=None, samples=None, ymin=-0.3, ymax=6.3,
    fontsize=25
):
    """
    Plot copy number profile from CovFrame[CopyNumber].

    Parameters
    ----------
    copy_number : str or pypgx.Archive
        Archive file or object with the semantic type CovFrame[CopyNumber].
    fitted : bool, default: False
        If True, show the fitted line as well.
    path : str, optional
        Create plots in this directory.
    samples : str or list, optional
        Specify which samples should be included for analysis by providing a
        text file (.txt, .tsv, .csv, or .list) containing one sample per
        line. Alternatively, you can provide a list of samples.
    ymin : float, default: -0.3
        Y-axis bottom.
    ymax : float, default: 6.3
        Y-axis top.
    fontsize : float, default: 25
        Text fontsize.
    """
    if isinstance(copy_number, str):
        copy_number = sdk.Archive.from_file(copy_number)

    copy_number.check_type('CovFrame[CopyNumber]')

    gene = copy_number.metadata['Gene']
    assembly = copy_number.metadata['Assembly']

    if samples is None:
        samples = copy_number.data.samples
    else:
        samples = common.parse_list_or_file(samples)
        copy_number = utils.filter_samples(copy_number, samples=samples)

    if fitted:
        processed_copy_number = utils._process_copy_number(copy_number)
    else:
        processed_copy_number = None

    with sns.axes_style('darkgrid'):
        for sample in samples:

            fig, [ax1, ax2] = plt.subplots(2, 1, figsize=(18, 12),
                gridspec_kw={'height_ratios': [1, 10]})

            ax1, ax2 = _plot_bam_copy_number_one(
                ax1, ax2, sample, copy_number, gene, assembly,
                processed_copy_number, ymin, ymax, fontsize
            )

            if path is None:
                output = f'{sample}.png'
            else:
                output = f'{path}/{sample}.png'

            plt.tight_layout()
            fig.savefig(output)
            plt.close()

def plot_bam_read_depth(
    read_depth, path=None, samples=None, ymin=None, ymax=None, fontsize=25
):
    """
    Plot copy number profile with BAM data.

    Parameters
    ----------
    read_depth : str or pypgx.Archive
        Archive file or object with the semantic type CovFrame[ReadDepth].
    path : str, optional
        Create plots in this directory.
    samples : str or list, optional
        Specify which samples should be included for analysis by providing a
        text file (.txt, .tsv, .csv, or .list) containing one sample per
        line. Alternatively, you can provide a list of samples.
    ymin : float, optional
        Y-axis bottom.
    ymax : float, optional
        Y-axis top.
    fontsize : float, default: 25
        Text fontsize.
    """

    if isinstance(read_depth, str):
        read_depth = sdk.Archive.from_file(read_depth)

    read_depth.check_type('CovFrame[ReadDepth]')

    if samples is None:
        samples = read_depth.data.samples
    else:
        samples = common.parse_list_or_file(samples)

    gene = read_depth.metadata['Gene']
    assembly = read_depth.metadata['Assembly']
    region = core.get_region(gene, assembly=assembly)
    chrom, start, end = common.parse_region(region)

    with sns.axes_style('darkgrid'):
        for sample in samples:

            fig, [ax1, ax2] = plt.subplots(2, 1, figsize=(18, 12), gridspec_kw={'height_ratios': [1, 10]})

            _plot_exons(gene, assembly, ax1)

            read_depth.data.plot_region(sample, ax=ax2, legend=False)

            ax2.set_xlim([start, end])
            ax2.set_ylim([ymin, ymax])
            ax2.set_xlabel('Coordinate (Mb)', fontsize=fontsize)
            ax2.set_ylabel('Read depth', fontsize=fontsize)
            ax2.tick_params(axis='both', which='major', labelsize=fontsize)
            ax2.ticklabel_format(axis='x', useOffset=False, scilimits=(6, 6))

            if path is None:
                output = f'{sample}.png'
            else:
                output = f'{path}/{sample}.png'

            plt.tight_layout()
            fig.savefig(output)
            plt.close()

def plot_cn_af(
    copy_number, imported_variants, path=None, samples=None, ymin=-0.3,
    ymax=6.3, fontsize=25
):
    """
    Plot both copy number profile and allele fraction profile in one figure.

    Parameters
    ----------
    copy_number : str or pypgx.Archive
        Archive file or object with the semantic type CovFrame[CopyNumber].
    imported_variants : str or pypgx.Archive
        Archive file or object with the semantic type VcfFrame[Imported] or
        VcfFrame[Consolidated].
    path : str, optional
        Create plots in this directory.
    samples : str or list, optional
        Specify which samples should be included for analysis by providing a
        text file (.txt, .tsv, .csv, or .list) containing one sample per
        line. Alternatively, you can provide a list of samples.
    ymin : float, default: -0.3
        Y-axis bottom.
    ymax : float, default: 6.3
        Y-axis top.
    fontsize : float, default: 25
        Text fontsize.
    """
    if isinstance(copy_number, str):
        copy_number = sdk.Archive.from_file(copy_number)

    copy_number.check_type('CovFrame[CopyNumber]')

    if isinstance(imported_variants, str):
        imported_variants = sdk.Archive.from_file(imported_variants)

    imported_variants.check_type(
        ['VcfFrame[Imported]', 'VcfFrame[Consolidated]'])

    sdk.compare_metadata('Gene', copy_number, imported_variants)
    sdk.compare_metadata('Assembly', copy_number, imported_variants)

    if samples is None:
        samples = copy_number.data.samples
    else:
        samples = common.parse_list_or_file(samples)
        copy_number = utils.filter_samples(copy_number, samples=samples)

    processed_copy_number = utils._process_copy_number(copy_number)

    gene = copy_number.metadata['Gene']
    assembly = copy_number.metadata['Assembly']

    with sns.axes_style('darkgrid'):
        for sample in samples:

            fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2, 2, figsize=(20, 10),
                gridspec_kw={'height_ratios': [1, 10]})

            ax1, ax3 = _plot_bam_copy_number_one(
                ax1, ax3, sample, copy_number, gene, assembly,
                processed_copy_number, ymin, ymax, fontsize
            )

            ax2, ax4 = _plot_vcf_allele_fraction_one(
                ax2, ax4, sample, imported_variants, gene, assembly, fontsize
            )

            if path is None:
                output = f'{sample}.png'
            else:
                output = f'{path}/{sample}.png'

            plt.tight_layout()
            fig.savefig(output)
            plt.close()

def plot_vcf_allele_fraction(
    imported_variants, path=None, samples=None, fontsize=25
):
    """
    Plot allele fraction profile with VCF data.

    Parameters
    ----------
    imported_variants : str or pypgx.Archive
        Archive file or object with the semantic type VcfFrame[Imported] or
        VcfFrame[Consolidated].
    path : str, optional
        Create plots in this directory.
    samples : str or list, optional
        Specify which samples should be included for analysis by providing a
        text file (.txt, .tsv, .csv, or .list) containing one sample per
        line. Alternatively, you can provide a list of samples.
    fontsize : float, default: 25
        Text fontsize.
    """
    if isinstance(imported_variants, str):
        imported_variants = sdk.Archive.from_file(imported_variants)

    imported_variants.check_type(
        ['VcfFrame[Imported]', 'VcfFrame[Consolidated]'])

    gene = imported_variants.metadata['Gene']
    assembly = imported_variants.metadata['Assembly']

    if samples is None:
        samples = imported_variants.data.samples
    else:
        samples = common.parse_list_or_file(samples)

    with sns.axes_style('darkgrid'):
        for sample in samples:

            fig, [ax1, ax2] = plt.subplots(2, 1, figsize=(18, 12), gridspec_kw={'height_ratios': [1, 10]})

            ax1, ax2 = _plot_vcf_allele_fraction_one(
                ax1, ax2, sample, imported_variants, gene, assembly, fontsize
            )

            if path is None:
                output = f'{sample}.png'
            else:
                output = f'{path}/{sample}.png'

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
    samples : str or list, optional
        Specify which samples should be included for analysis by providing a
        text file (.txt, .tsv, .csv, or .list) containing one sample per
        line. Alternatively, you can provide a list of samples.
    ymin : float, optional
        Y-axis bottom.
    ymax : float, optional
        Y-axis top.
    """
    vf = pyvcf.VcfFrame.from_file(vcf)

    region = core.get_region(gene, assembly=assembly)
    chrom, start, end = common.parse_region(region)

    if samples is None:
        samples = cf.samples
    else:
        samples = common.parse_list_or_file(samples)

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
