Changelog
*********

0.9.0 (in development)
----------------------

* Add 1KGP reference haplotype panels for GRCh37 for the 17 recently added genes.
* Add 1KGP reference haplotype panels for GRCh38 for all target genes.
* Add GRCh37 CNV callers for UGT1A4.
* Add GRCh38 CNV callers for all target genes.
* Update commands :command:`run-ngs-pipeline` and :command:`run-chip-pipeline` to support GRCh38.
* Update the **GeT-RM WGS** tutorial to include more genes and also GRCh38.
* Make the profiles (e.g. copy number) look prettier.
* Rename :meth:`sdk.utils.Archive.check` method to :meth:`sdk.utils.Archive.check_type`.
* Add new method :meth:`sdk.utils.Archive.check_metadata`.
* Add new error ``sdk.utils.IncorrectMetadataError``.
* Update :command:`run-ngs-pipeline` command to check input files more vigorously.
* Add new method :meth:`sdk.utils.compare_metadata`.

0.8.0 (2021-11-20)
------------------

* Update :meth:`api.core.sort_alleles` method to also sort alleles by name for genes that do not use the star allele nomenclature (e.g. the DPYD gene).
* Add new method :meth:`api.core.is_legit_allele`.
* Update :meth:`api.core.predict_phenotype` method to first check if the two alleles are legit.
* Add new genes: ABCB1, CYP1A1, CYP1B1, CYP4A11, CYP4A22, CYP4B1, CYP17A1, CYP19A1, G6PD, IFNL3, POR, PTGIS, SLCO1B3, SULT1A1, TBXAS1, UGT1A4, XPC.

0.7.0 (2021-10-23)
------------------

* Fix minor bug in :meth:`api.core.predict_phenotype` when specified diplotype is not present in diplotype table.
* Dissolve **Database of Pharmacogenomic Structural Variants (DPSV)** page and move its SV data to **Genes** page.
* Add new method :meth:`api.core.get_variant_impact`.
* Update :meth:`api.utils.sort_alleles` method to give priority to alleles that impact protein coding when breaking ties (i.e. alleles have the same functional status and same number of variants).
* Update CNV caller for SLC22A2 and UGT2B15 genes.
* Rename ``--chr-prefix`` argument in :command:`create-regions-bed` to ``--add-chr-prefix``.
* Add ``--samples`` argument to :command:`run-ngs-pipeline` command.
* Add new command :command:`compare-genotypes`.
* Update :meth:`api.genotype.call_genotypes` method to assume the samples have no SV when CNV calls are not provided even if the target gene is known to have SV.
* Add new command :command:`run-chip-pipeline`.
* Fix minor bug in :command:`estimate-phase-beagle` command on not properly exiting the program even though there was an error raised by Beagle.
* Update :meth:`api.utils.create_consolidated_vcf` method to check synonymous variants as well when performing phase-extension algorithm.
* Update :command:`run-ngs-pipeline` command to output a warning when user provides CovFrame[DepthOfCoverage] even though target gene does not have any star alleles defined by SVs.
* Add new argument ``--fontsize`` argument to :command:`plot-bam-copy-number` command.
* Remove ``--ymin`` and ``--ymax`` arguments from :command:`plot-vcf-allele-fraction` command.
* Update ``--ymin`` and ``--ymax`` arguments of :command:`plot-bam-copy-number` command to have a default value.
* Add new command :command:`plot-cn-af`.
* Update :command:`run-ngs-pipeline` command to output a warning when user provides a VCF file even though target gene does not have any star alleles defined by SNVs/indels.
* Update aesthetics of copy number profile and allele fraction profile.
* Add new method :meth:`api.utils.count_alleles`.
* Update variant information for following alleles: CYP2A6\*35, UGT1A1\*28, UGT1A1\*37.

0.6.0 (2021-10-09)
------------------

* :issue:`25`: Add new extension ``sphinx-issues`` to Read the Docs.
* :issue:`26`: Add new extension ``sphinx.ext.linkcode`` to Read the Docs.
* Add ``by`` argument to :meth:`api.utils.sort_alleles` method. When ``by='name'`` it will sort star alleles by allele number.
* Update :command:`call-genotypes` command to output genotypes with number-sorted alleles (e.g. '\*4/\*10' instead of '\*10/\*4').
* Add new semantic type ``SampleTable[Phenotypes]``.
* Add new method :meth:`api.utils.call_phenotypes`.
* Add new command :command:`call-phenotypes`.
* Add ``--phenotypes`` argument  to :command:`combine-results` command.
* Deprecate :meth:`api.utils.load_control_table` method.
* Split ``api.utils`` submodule into two submodules ``api.utils`` and ``api.core``.
* Update :command:`run-ngs-pipeline` command to include phenotype calling step.
* Update :command:`plot-bam-copy-number` command to run faster when ``--samples`` argument is used.
* Change 'Unassigned' genotype to 'Indeterminate' genotype.
* Add new method :meth:`api.core.get_variant_synonyms`.
* Update :meth:`api.core.list_variants` method to accept multiple star alleles.
* Update :command:`predict-alleles` command to support multiallelic variants.
* Update :meth:`api.utils.sort_alleles` method to give priority to non-reference or non-default alleles when breaking ties (i.e. alleles have the same functional status and same number of variants).
* Update variant information for following alleles: CYP2D6\*122, CYP2D6\*127, CYP2D6\*139.

0.5.0 (2021-10-02)
------------------

* Update :command:`create-read-depth-tsv` command to automatically detect ``chr`` string in input BAM.
* Add ``sdk.utils.parse_input_bams`` method.
* Add the 1000 Genomes Project reference haplotype panel for GRCh37. When estimating haplotype phase of observed variants, users are no longer needed to download and specify a panel. GRCh38 support will follow in a future release.
* Rename command :command:`create-read-depth-tsv` to :command:`prepare-depth-of-coverage`.
* Add ``bed`` argument to :command:`prepare-depth-of-coverage` command.
* Update :command:`prepare-depth-of-coverage` command to output archive file instead of TSV file.
* Update :command:`import-read-depth` command to accept archive file as input instead of TSV file.
* Add ``fitted`` argument to :command:`plot-bam-copy-number` command.
* From now on, missing copy number will be imputed with forward filling instead of column median.
* Update :command:`predict-cnv` command to support a user-defined CNV caller.
* Add **Database of Pharmacogenomic Structural Variants (DPSV)** page.
* Update :command:`predict-alleles` command to output variant data even for alleles in ``AlternativePhase`` column.
* Update :command:`create-consolidated-vcf` command to mark phased variants with 'Phased' in ``INFO`` column in VCF.
* Update the allele table.
* Update :meth:`api.utils.list_alleles` method to be able to only list alleles carrying specified variant(s) as a part of definition.
* Add ``mode`` argument to :meth:`api.utils.list_variants` method.
* Update :command:`create-consolidated-vcf` command to implement phase-extension algorithm.
* Remove ``SO`` and ``Type`` columns from the variant table.
* Update :class:`api.genotype.GSTM1Genotyper` class.
* Add ``NotTargetGeneError`` error.
* Add new method ``api.utils.is_target_gene``.
* Update :command:`run-ngs-pipeline` command to check whether input gene is one of the target genes before attempting to run the pipeline.
* Update variant information for following alleles: CYP1A2\*1C, CYP1A2\*1F, CYP1A2\*1K, CYP1A2\*1L, CYP2B6\*17, CYP2D6\*15, CYP2D6\*21, SLCO1B1\*S1, SLCO1B1\*S2.

0.4.1 (2021-09-21)
------------------

* Initial release.
