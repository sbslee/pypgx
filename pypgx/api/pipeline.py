"""
The pipeline submodule is used to provide convenient methods that combine
multiple PyPGx actions and automatically handle semantic types.
"""

import shutil
import os

from .. import sdk

from . import utils, plot, genotype

def run_ngs_pipeline(
    gene, output, variants=None, depth_of_coverage=None,
    control_statistics=None, panel=None,
    force=False, do_not_plot_copy_number=False,
    do_not_plot_allele_fraction=False
):
    """
    Run WGS pipeline for a target gene.

    Parameters
    ----------
    gene : str
        Target gene.
    output : str
        Output directory.
    variants : str, optional
        VCF file (zipped or unzipped).
    depth_of_coverage : str, optional
        Archive file or object with the semantic type
        CovFrame[DepthOfCoverage].
    control_statistics : str or pypgx.Archive, optional
        Archive file or object with the semantic type SampleTable[Statistics].
    panel : str, optional
        Reference haplotype panel.
    force : bool, default : False
        Overwrite output directory if it already exists.
    do_not_plot_copy_number : bool, default: False
        Do not plot copy number profile.
    do_not_plot_allele_fraction : bool, default: False
        Do not plot allele fraction profile.
    """
    if not utils.is_target_gene(gene):
        raise utils.NotTargetGeneError(gene)

    if os.path.exists(output) and force:
        shutil.rmtree(output)

    os.mkdir(output)

    gene_table = utils.load_gene_table()

    alleles = None
    cnv_calls = None

    if gene_table[gene_table.Gene == gene].Variants.values[0] and variants is not None:
        imported_variants = utils.import_variants(gene, variants)
        imported_variants.to_file(f'{output}/imported-variants.zip')
        phased_variants = utils.estimate_phase_beagle(imported_variants, panel=panel)
        phased_variants.to_file(f'{output}/phased-variants.zip')
        consolidated_variants = utils.create_consolidated_vcf(imported_variants, phased_variants)
        consolidated_variants.to_file(f'{output}/consolidated-variants.zip')
        alleles = utils.predict_alleles(consolidated_variants)
        alleles.to_file(f'{output}/alleles.zip')
        if not do_not_plot_allele_fraction:
            os.mkdir(f'{output}/allele-fraction-profile')
            plot.plot_vcf_allele_fraction(
                imported_variants, path=f'{output}/allele-fraction-profile'
            )

    if gene_table[gene_table.Gene == gene].SV.values[0] and depth_of_coverage is not None:
        if isinstance(depth_of_coverage, str):
            depth_of_coverage = sdk.Archive.from_file(depth_of_coverage)

        depth_of_coverage.check('CovFrame[DepthOfCoverage]')

        if control_statistics is None:
            raise ValueError('SV detection requires SampleTable[Statistcs]')

        if isinstance(control_statistics, str):
            control_statistics = sdk.Archive.from_file(control_statistics)

        control_statistics.check('SampleTable[Statistics]')

        if depth_of_coverage.metadata['Platform'] != control_statistics.metadata['Platform']:
            raise ValueError('Different platforms detected')

        read_depth = utils.import_read_depth(gene, depth_of_coverage)
        read_depth.to_file(f'{output}/read-depth.zip')
        copy_number = utils.compute_copy_number(read_depth, control_statistics)
        copy_number.to_file(f'{output}/copy-number.zip')
        cnv_calls = utils.predict_cnv(copy_number)
        cnv_calls.to_file(f'{output}/cnv-calls.zip')
        if not do_not_plot_copy_number:
            os.mkdir(f'{output}/copy-number-profile')
            plot.plot_bam_copy_number(
                copy_number, path=f'{output}/copy-number-profile', ymin=0,
                ymax=6
            )

    genotypes = genotype.call_genotypes(alleles=alleles, cnv_calls=cnv_calls)
    genotypes.to_file(f'{output}/genotypes.zip')
    results = utils.combine_results(genotypes=genotypes, alleles=alleles, cnv_calls=cnv_calls)
    results.to_file(f'{output}/results.zip')
