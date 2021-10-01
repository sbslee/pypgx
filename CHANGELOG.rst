Changelog
*********

0.5.0 (in development)
----------------------

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
* Update variant information for following alleles: CYP2B6*17, CYP2D6*15, CYP2D6*21, SLCO1B1*S1, SLCO1B1*S2.

0.4.1 (2021-09-21)
------------------

* Initial release.
