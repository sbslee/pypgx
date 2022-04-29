"""
The pipeline submodule is used to provide convenient methods that combine
multiple PyPGx actions and automatically handle semantic types.
"""

import shutil
import os
import warnings

from .. import sdk

from . import utils, plot, genotype, core

def run_chip_pipeline(
    gene, output, variants, assembly='GRCh37', panel=None, impute=False,
    force=False, samples=None, exclude=False
):
    """
    Run genotyping pipeline for chip data.

    Parameters
    ----------
    gene : str
        Target gene.
    output : str
        Output directory.
    variants : str
        Input VCF file must be already BGZF compressed (.gz) and indexed
        (.tbi) to allow random access. Statistical haplotype phasing will be
        skipped if input VCF is already fully phased.
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.
    panel : str, optional
        VCF file corresponding to a reference haplotype panel (compressed or
        uncompressed). By default, the 1KGP panel in the ``~/pypgx-bundle``
        directory will be used.
    impute : bool, default: False
        If True, perform imputation of missing genotypes.
    force : bool, default : False
        Overwrite output directory if it already exists.
    samples : str or list, optional
        Subset the VCF for specified samples. This can be a text file (.txt,
        .tsv, .csv, or .list) containing one sample per line. Alternatively,
        you can provide a list of samples.
    exclude : bool, default: False
        If True, exclude specified samples.
    """
    if not core.is_target_gene(gene):
        raise sdk.utils.NotTargetGeneError(gene)

    if os.path.exists(output) and force:
        shutil.rmtree(output)

    os.mkdir(output)

    imported_variants = utils.import_variants(gene, variants,
        assembly=assembly, platform='Chip', samples=samples, exclude=exclude)
    imported_variants.to_file(f'{output}/imported-variants.zip')

    # Skip statistical phasing if input VCF is already fully phased.
    if imported_variants.type == 'VcfFrame[Consolidated]':
        consolidated_variants = imported_variants
    else:
        phased_variants = utils.estimate_phase_beagle(
            imported_variants, panel=panel)
        phased_variants.to_file(f'{output}/phased-variants.zip')
        consolidated_variants = utils.create_consolidated_vcf(
            imported_variants, phased_variants)
        consolidated_variants.to_file(
            f'{output}/consolidated-variants.zip')

    alleles = utils.predict_alleles(consolidated_variants)
    alleles.to_file(f'{output}/alleles.zip')
    genotypes = genotype.call_genotypes(alleles=alleles)
    genotypes.to_file(f'{output}/genotypes.zip')
    phenotypes = utils.call_phenotypes(genotypes)
    phenotypes.to_file(f'{output}/phenotypes.zip')
    results = utils.combine_results(
        genotypes=genotypes, phenotypes=phenotypes, alleles=alleles
    )
    results.to_file(f'{output}/results.zip')

def run_long_read_pipeline(
    gene, output, variants, assembly='GRCh37', force=False, samples=None,
    exclude=False
):
    """
    Run genotyping pipeline for long-read sequencing data.

    Parameters
    ----------
    gene : str
        Target gene.
    output : str
        Output directory.
    variants : str
        Input VCF file must be already BGZF compressed (.gz) and indexed
        (.tbi) to allow random access.
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.
    force : bool, default : False
        Overwrite output directory if it already exists.
    samples : str or list, optional
        Subset the VCF for specified samples. This can be a text file (.txt,
        .tsv, .csv, or .list) containing one sample per line. Alternatively,
        you can provide a list of samples.
    exclude : bool, default: False
        If True, exclude specified samples.
    """
    if not core.is_target_gene(gene):
        raise sdk.utils.NotTargetGeneError(gene)

    if os.path.exists(output) and force:
        shutil.rmtree(output)

    os.mkdir(output)

    consolidated_variants = utils.import_variants(gene, variants,
        assembly=assembly, platform='LongRead', samples=samples,
        exclude=exclude)
    consolidated_variants.to_file(f'{output}/consolidated-variants.zip')
    alleles = utils.predict_alleles(consolidated_variants)
    alleles.to_file(f'{output}/alleles.zip')
    genotypes = genotype.call_genotypes(alleles=alleles)
    genotypes.to_file(f'{output}/genotypes.zip')
    phenotypes = utils.call_phenotypes(genotypes)
    phenotypes.to_file(f'{output}/phenotypes.zip')
    results = utils.combine_results(
        genotypes=genotypes, phenotypes=phenotypes, alleles=alleles
    )
    results.to_file(f'{output}/results.zip')

def run_ngs_pipeline(
    gene, output, variants=None, depth_of_coverage=None,
    control_statistics=None, platform='WGS', assembly='GRCh37', panel=None,
    force=False, samples=None, exclude=False, samples_without_sv=None,
    do_not_plot_copy_number=False, do_not_plot_allele_fraction=False,
    cnv_caller=None
):
    """
    Run genotyping pipeline for NGS data.

    During copy number analysis, if the input data is targeted sequencing,
    the method will apply inter-sample normalization using summary statistics
    across all samples. For best results, it is recommended to specify known
    samples without SV using ``samples_without_sv``.

    Parameters
    ----------
    gene : str
        Target gene.
    output : str
        Output directory.
    variants : str, optional
        Input VCF file must be already BGZF compressed (.gz) and indexed
        (.tbi) to allow random access. Statistical haplotype phasing will be
        skipped if input VCF is already fully phased.
    depth_of_coverage : str, optional
        Archive file or object with the semantic type
        CovFrame[DepthOfCoverage].
    control_statistics : str or pypgx.Archive, optional
        Archive file or object with the semantic type SampleTable[Statistics].
    platform : {'WGS', 'Targeted'}, default: 'WGS'
        Genotyping platform.
    assembly : {'GRCh37', 'GRCh38'}, default: 'GRCh37'
        Reference genome assembly.
    panel : str, optional
        VCF file corresponding to a reference haplotype panel (compressed or
        uncompressed). By default, the 1KGP panel in the ``~/pypgx-bundle``
        directory will be used.
    force : bool, default : False
        Overwrite output directory if it already exists.
    samples : str or list, optional
        Subset the VCF for specified samples. This can be a text file (.txt,
        .tsv, .csv, or .list) containing one sample per line. Alternatively,
        you can provide a list of samples.
    exclude : bool, default: False
        If True, exclude specified samples.
    samples_without_sv : list, optional
        List of known samples without SV.
    do_not_plot_copy_number : bool, default: False
        Do not plot copy number profile.
    do_not_plot_allele_fraction : bool, default: False
        Do not plot allele fraction profile.
    cnv_caller : str or pypgx.Archive, optional
        Archive file or object with the semantic type Model[CNV]. By default,
        a pre-trained CNV caller in the ``~/pypgx-bundle`` directory will be
        used.
    """
    if not core.is_target_gene(gene):
        raise sdk.utils.NotTargetGeneError(gene)

    gene_table = core.load_gene_table()
    small_var = gene_table[gene_table.Gene == gene].Variants.values[0]
    large_var = gene_table[gene_table.Gene == gene].SV.values[0]

    if not small_var and variants is not None:
        message = (
            'User provided a VCF file even though the target gene does '
            'not have any star alleles defined by SNVs/indels. PyPGx will '
            'ignore it.'
        )
        warnings.warn(message)

    if not large_var and depth_of_coverage is not None:
        message = (
            'User provided CovFrame[DepthOfCoverage] even though the '
            'target gene does not have any star alleles defined by SVs. '
            'PyPGx will ignore it.'
        )
        warnings.warn(message)

    if not large_var and control_statistics is not None:
        message = (
            'User provided SampleTable[Statistics] even though the '
            'target gene does not have any star alleles defined by SVs. '
            'PyPGx will ignore it.'
        )
        warnings.warn(message)

    alleles = None
    cnv_calls = None

    if os.path.exists(output) and force:
        shutil.rmtree(output)

    os.mkdir(output)

    if small_var and variants is not None:
        imported_variants = utils.import_variants(gene, variants,
            assembly=assembly, platform=platform, samples=samples,
            exclude=exclude)
        imported_variants.to_file(f'{output}/imported-variants.zip')

        # Skip statistical phasing if input VCF is already fully phased.
        if imported_variants.type == 'VcfFrame[Consolidated]':
            consolidated_variants = imported_variants
        else:
            phased_variants = utils.estimate_phase_beagle(
                imported_variants, panel=panel)
            phased_variants.to_file(f'{output}/phased-variants.zip')
            consolidated_variants = utils.create_consolidated_vcf(
                imported_variants, phased_variants)
            consolidated_variants.to_file(
                f'{output}/consolidated-variants.zip')

        alleles = utils.predict_alleles(consolidated_variants)
        alleles.to_file(f'{output}/alleles.zip')

        if not do_not_plot_allele_fraction:
            os.mkdir(f'{output}/allele-fraction-profile')
            plot.plot_vcf_allele_fraction(
                imported_variants, path=f'{output}/allele-fraction-profile'
            )

    if large_var and depth_of_coverage is not None:
        if isinstance(depth_of_coverage, str):
            depth_of_coverage = sdk.Archive.from_file(depth_of_coverage)

        depth_of_coverage.check_type('CovFrame[DepthOfCoverage]')
        depth_of_coverage.check_metadata('Platform', platform)
        depth_of_coverage.check_metadata('Assembly', assembly)

        if control_statistics is None:
            raise ValueError('SV detection requires SampleTable[Statistics]')

        if isinstance(control_statistics, str):
            control_statistics = sdk.Archive.from_file(control_statistics)

        if samples is not None:
            control_statistics = utils.filter_samples(control_statistics,
                samples=samples, exclude=exclude)

        control_statistics.check_type('SampleTable[Statistics]')
        control_statistics.check_metadata('Platform', platform)
        control_statistics.check_metadata('Assembly', assembly)

        read_depth = utils.import_read_depth(gene, depth_of_coverage,
            samples=samples, exclude=exclude)
        read_depth.to_file(f'{output}/read-depth.zip')
        copy_number = utils.compute_copy_number(read_depth,
            control_statistics, samples_without_sv=samples_without_sv)
        copy_number.to_file(f'{output}/copy-number.zip')
        cnv_calls = utils.predict_cnv(copy_number, cnv_caller=cnv_caller)
        cnv_calls.to_file(f'{output}/cnv-calls.zip')
        if not do_not_plot_copy_number:
            os.mkdir(f'{output}/copy-number-profile')
            plot.plot_bam_copy_number(
                copy_number, path=f'{output}/copy-number-profile'
            )

    genotypes = genotype.call_genotypes(alleles=alleles, cnv_calls=cnv_calls)
    genotypes.to_file(f'{output}/genotypes.zip')
    phenotypes = utils.call_phenotypes(genotypes)
    phenotypes.to_file(f'{output}/phenotypes.zip')
    results = utils.combine_results(
        genotypes=genotypes, phenotypes=phenotypes, alleles=alleles,
        cnv_calls=cnv_calls
    )
    results.to_file(f'{output}/results.zip')
