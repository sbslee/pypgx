from .api.utils import (
    build_definition_table,
    create_consolidated_vcf,
    get_default_allele,
    get_function,
    get_paralog,
    get_priority,
    get_region,
    get_score,
    has_definition,
    has_phenotype,
    has_score,
    list_alleles,
    list_functions,
    list_genes,
    list_phenotypes,
    list_variants,
    load_allele_table,
    load_control_table,
    load_diplotype_table,
    load_equation_table,
    load_gene_table,
    load_phenotype_table,
    load_variant_table,
    predict_alleles,
    predict_phenotype,
    predict_score,
    sort_alleles,
)

from .api.plot import (
    plot_allele_fraction,
    plot_copy_number,
    plot_read_depth,
)

del api
