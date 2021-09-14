import shutil
import os

from . import utils, plot, genotype

def run_ngs_pipeline(
    gene, output, vcf=None, panel=None, tsv=None, control=None, force=False,
    plot_copy_number=True
):
    """
    Run NGS pipeline for target gene.
    """
    if os.path.exists(output) and force:
        shutil.rmtree(output)

    os.mkdir(output)

    gene_table = utils.load_gene_table()

    alleles = None
    cnv_calls = None

    if gene_table[gene_table.Gene == gene].Variants.values[0] and vcf is not None:
        imported_variants = utils.import_variants(gene, vcf)
        imported_variants.to_file(f'{output}/imported-variants.zip')
        phased_variants = utils.estimate_phase_beagle(imported_variants, panel)
        phased_variants.to_file(f'{output}/phased-variants.zip')
        consolidated_variants = utils.create_consolidated_vcf(imported_variants, phased_variants)
        consolidated_variants.to_file(f'{output}/consolidated-variants.zip')
        alleles = utils.predict_alleles(consolidated_variants)
        alleles.to_file(f'{output}/alleles.zip')

    if gene_table[gene_table.Gene == gene].SV.values[0] and tsv is not None:
        if control is None:
            raise ValueError('CovFrame[ReadDepth] requires SampleTable[Statistcs]')
        read_depth = utils.import_read_depth(gene, tsv)
        read_depth.to_file(f'{output}/read-depth.zip')
        copy_number = utils.compute_copy_number(read_depth, control)
        copy_number.to_file(f'{output}/copy-number.zip')
        cnv_calls = utils.predict_cnv(copy_number)
        cnv_calls.to_file(f'{output}/cnv-calls.zip')
        if plot_copy_number:
            os.mkdir(f'{output}/plots')
            plot.plot_bam_copy_number(
                copy_number, path=f'{output}/plots', ymin=0, ymax=6
            )

    genotypes = genotype.call_genotypes(alleles=alleles, cnv_calls=cnv_calls)

    genotypes.to_file(f'{output}/genotypes.zip')

    results = utils.combine_results(genotypes=genotypes, alleles=alleles, cnv_calls=cnv_calls)

    results.to_file(f'{output}/results.zip')
