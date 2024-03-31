Changelog
*********

0.25.0 (in development)
-----------------------

0.24.0 (2024-03-31)
-------------------

* Change G6PD haplotype name from ``*MALE`` to ``MALE``, as the gene does not utilize star allele nomenclature.
* Update G6PD haplotype nomenclature in accordance with the latest CPIC version (2024-01-07).
* :issue:`119`: Add G6PD phenotype data in accordance with the latest CPIC version (2024-01-07).
* Update CYP4F2 haplotype nomenclature in accordance with the latest PharmVar version (2024-01-27).

0.23.0 (2023-12-24)
-------------------

* :issue:`115`: Add new genes COMT and MTHFR (thanks `@nbiesot <https://github.com/nbiesot>`__).
* :issue:`117`: Fix major bug in :meth:`api.utils.estimate_phase_beagle` method. The bug was introduced in version 0.22.0 (:issue:`113`) and carries significant implications, especially if the user's input VCF contains sample names that overlap with the reference panel or if there are differences in chromosome annotation compared to the reference panel (e.g. 'chr22' vs. '22').

0.22.0 (2023-12-11)
-------------------

* :issue:`100`: Add new method :meth:`sdk.utils.get_bundle_path` to enable customization of the ``pypgx-bundle`` directory's location instead of the user's home directory.
* :issue:`114`: Fix bug in :meth:`api.core.get_recommendation` method where string ``'None'`` was treated as missing value by ``pandas.read_csv`` version 2.0 or higher.
* :issue:`113`: Fix bug in :meth:`api.utils.estimate_phase_beagle` method where Beagle's expectation-maximization algorithm estimated a parameter value that was outside the permitted range.

0.21.0 (2023-08-25)
-------------------

* :issue:`96`: Fix bug in :command:`compute-copy-number` command where it threw an error "Different sample sets found" when the sample name consists of only numbers with one or more leading zeros.
* :issue:`97`: Fix bug in IFNL3 genotyping where certain variants were not considered because they were located outside of the target region.
* :issue:`99`: Refresh the CPIC table for :meth:`api.core.load_cpic_table` method. The table's latest update was on February 6, 2022.
* :issue:`63`: Fix bug in :meth:`api.utils.estimate_phase_beagle` method where Beagle throws an error "Window has only one position" even when multiple overlapping variants exist.

0.20.0 (2023-01-12)
-------------------

* :issue:`73`: Fix bug in :command:`run-ngs-pipeline` command where empty VCF (i.e. no variants were found in the target gene) was causing error when plotting allele fraction. From now on, a warning will be produced telling users to specify ``--do-not-plot-allele-fraction`` to suppress this warning.
* :issue:`78`: Fix bug in :command:`compute-copy-number` command where it threw an error "Different sample sets found" when the sample name consists of only numbers.

0.19.0 (2022-09-13)
-------------------

* Add new method :meth:`api.core.has_sv`.
* Update :meth:`api.core.sort_alleles` method to handle ``'Indeterminate'`` haplotype call when ``by='name'``.
* Update :meth:`api.utils.estimate_phase_beagle` method to handle situations where there are overlapping samples between input VCF and reference panel -- i.e. users are no longer required to change sample names. Before this update, the Beagle program would throw an error if there were overlapping samples (e.g. 1KGP samples). From now on, offending samples will be temporarily renamed before statistical phasing.
* Add new methods :meth:`api.core.load_recommendation_table` and :meth:`api.core.get_recommendation`.

0.18.0 (2022-08-12)
-------------------

* PyPGx now has a citation! Please refer to the publication "`ClinPharmSeq: A targeted sequencing panel for clinical pharmacogenetics implementation <https://doi.org/10.1371/journal.pone.0272129>`__" by Lee et al., 2022 (Steven is the first author). Fore more details, see the Citation section in README.
* Update phenotype data and star allele nomenclature for SLCO1B1 in accordance with the latest PharmVar version (v5.2.1). Note that SLCO1B1 was first formally added to PharmVar v5.1 on October 14, 2021. For more details, please refer to the publication "`PharmVar GeneFocus: SLCO1B1 <https://doi.org/10.1002/cpt.2705>`__" by Ramsey et al., 2022 (Steven is a co-author of this paper by the way) and the change log in `the PharmVar SLCO1B1 page <https://www.pharmvar.org/gene/SLCO1B1>`__. The PharmVar-developed SLCO1B1 nomenclature has been incorporated by CPIC 2022 guideline on statin-associated musculoskeletal symptoms.
* Remove duplicate DYPD entry from ``phenotype-table.csv`` file (i.e. Poor Metabolizer).
* Fix major bug in :command:`run-chip-pipeline` command where ``--impute`` argument is essentially ignored.
* :issue:`68`: Fix bug in :meth:`api.utils.estimate_phase_beagle` method when there are no overlapping variants between input VCF and reference panel.
* :issue:`68`: Update :meth:`api.utils.estimate_phase_beagle` method to warn when statistical phasing is skipped.
* :issue:`68`: Upgrade Beagle version from v5.2 (beagle.28Jun21.220.jar) to v5.4 (beagle.22Jul22.46e.jar) due to a bug in v5.2.
* :issue:`68`: Update :meth:`api.utils.estimate_phase_beagle` method to filter out variants with improper allele ('I', 'D', 'N', '.'). Note that this issue is specific to chip data.
* :issue:`68`: Update :meth:`api.utils.import_variants` method to handle input VCF with duplicate variants. Basically, it will warn the user about it and and only keep the first record. This issue seems to occur frequently with chip data.

0.17.0 (2022-07-12)
-------------------

* :issue:`63`: Fix bug in :meth:`api.utils.estimate_phase_beagle` method when there is only one variant in input VCF and Beagle throws an error.
* Update :command:`compare-genotypes` command to print the entire discordant calls when ``--verbose`` is used.
* Update :command:`compute-copy-number` command to ensure that the samples in CovFrame[ReadDepth] and SampleTable[Statistics] are in the same order.
* :issue:`64`: Update :meth:`api.utils.import_variants` method to 'diploidize' the input VCF when the target gene is G6PD. This is because some variant callers output haploid genotypes for males for the X chromosome, interfering with downstream analyses.
* Remove unnecessary optional argument ``assembly`` from :meth:`api.core.get_ref_allele`.

0.16.0 (2022-06-08)
-------------------

* Add new optional argument ``--comparison-table`` to :command:`train-cnv-caller` and :command:`test-cnv-caller` commands.
* Update :meth:`sdk.utils.add_cn_samples` and :meth:`sdk.utils.simulate_copy_number` methods to check input files more rigorously.
* Update :meth:`api.utils.test_cnv_caller` and :meth:`api.utils.train_cnv_caller` methods to accept the latest format of SampleTable[CNVCalls] as input.
* Update plotting methods to optionally return a list of :class:`matplotlib.figure.Figure` objects for API users (e.g. Jupyter Notebook): :meth:`api.plot.plot_bam_copy_number`, :meth:`api.plot.plot_bam_read_depth`, :meth:`api.plot.plot_cn_af`, :meth:`api.plot.plot_vcf_allele_fraction`, :meth:`api.plot.plot_vcf_read_depth`.
* :issue:`61`: Fix bug in commands :command:`compute-control-statistics`, :command:`compute-target-depth`, and :command:`prepare-depth-of-coverage` when a BED file is provided by user.
* Improve CNV caller for CYP2A6, CYP2B6, CYP2D6, CYP2E1, GSTM1, SLC22A2, SULT1A1, UGT1A4, UGT2B15, UGT2B17.
* Add new CNV call for CYP2A6: ``Unknown1``, ``Hybrid7``, ``Tandem2``.
* Add new CNV calls for CYP2B6: ``Tandem1``, ``PartialDup1``, ``PartialDup2``, ``ParalogWholeDel1``.
* Add new CNV call for CYP2D6: ``WholeDel1+Tandem3``. Also, remove ``PseudogeneDownstreamDel``.
* Add new CNV calls for CYP2E1: ``WholeDel1`` and ``WholeDup1+PartialDup1``.
* Add new CNV call for SLC22A2: ``NoncodingDel1Hom``.
* Add new CNV call for SULT1A1: ``Unknown2``, ``Unknown3``, ``Unknown4``.
* Add new CNV call for UGT1A4: ``NoncodingDel1Hom``.
* Add new CNV call for UGT2B15: ``PartialDup2``.
* Add new CNV call for UGT2B17: ``PartialDel2``. Also, define a new star allele ``*S3`` for ``PartialDel3``.
* :issue:`59`: Update CNV labels.

0.15.0 (2022-05-03)
-------------------

* Add new optional arguments ``--genes`` and ``--exclude`` to :command:`prepare-depth-of-coverage` command.
* Add new command :command:`slice-bam`.
* Add new command :command:`print-data`.
* Fix typo "statistcs" to "statistics" throughout the package.
* Update :meth:`sdk.utils.simulate_copy_number` method to automatically handle duplicate sample names.
* Improve CNV caller for CYP2A6, CYP2B6, CYP2D6, CYP2E1, GSTM1, SLC22A2, SULT1A1, UGT1A4, UGT2B15, UGT2B17.
* Add new CNV calls for CYP2A6: ``Deletion2Hom``, ``Hybrid5``, ``Hybrid6``, ``PseudogeneDeletion``.
* Add new CNV call for CYP2D6: ``Tandem2F``.
* Add new CNV call for GSTM1: ``Normal,Deletion2``.
* Add new CNV call for SULT1A1: ``Unknown1``.
* Add new CNV call for UGT2B17: ``Deletion,PartialDeletion3``.

0.14.0 (2022-04-03)
-------------------

* :issue:`49`: Add new gene ABCG2.
* :issue:`50`: Improve algorithm for whole gene duplication detection. This minor update will handle some rare edge cases.
* :issue:`53`: Update CYP2B6\*17 variants to have synonyms. Update :meth:`api.core.get_variant_synonyms` and :meth:`api.utils.predict_alleles` methods to allow mapping of single variant to multiple synonyms.
* :issue:`54`: Add new command :command:`create-input-vcf`.
* Fix minor error in ``gene-table.csv``. Two control genes EGFR and VDR incorrectly had ``TRUE`` for ``Variants`` column. They were changed to ``FALSE``.
* Change the three paralogs in ``gene-table.csv`` (CYP2A7, CYP2B7, and CYP2D7) to have ``FALSE`` for ``SV`` column.
* Add new optional arguments to :command:`create-regions-bed` command: ``--target-genes``, ``--var-genes``, ``--genes``, ``--exclude``.
* Improve CNV caller for CYP2A6, CYP2B6, CYP2D6, CYP2E1, GSTM1, SLC22A2, SULT1A1, UGT1A4, UGT2B15, UGT2B17.
* Add new CNV calls for CYP2A6: ``Hybrid2Hom``, ``Hybrid4``, ``Tandem``.
* Add new CNV calls for CYP2D6: ``Tandem4``, ``PseudogeneDownstreamDel``, ``Unknown2``.
* Add a new CNV call for CYP2E1: ``Multiplication2``.
* Add new CNV calls for GSTM1: ``PartialDuplication`` and ``DeletionHet,Deletion2``.
* Add a new CNV call for SLC22A2: ``PartialDuplication``.
* Add a new CNV call for SULT1A1: ``DeletionHom``.
* Add new CNV calls for UGT2B15: ``Deletion2``, ``Duplication``, ``PartialDuplication``.
* Add a new CNV call for UGT2B17: ``Normal,PartialDeletion3``.

0.13.0 (2022-03-01)
-------------------

* Add new genotyping platform, ``LongRead``, to :command:`import-variants` command.
* Add new command :command:`run-long-read-pipeline`.
* Remove ``Code`` column from ``cnv-table.csv`` file. From now on, CNV codes will be generated on the fly.
* Add new method :meth:`api.core.load_cpic_table`.
* Move following errors from ``api.core`` submodule to ``sdk.utils`` submodule: :class:`AlleleNotFoundError`, :class:`GeneNotFoundError`, :class:`NotTargetGeneError`, :class:`PhenotypeNotFoundError`, :class:`VariantNotFoundError`.
* Combine optional arguments ``--bam`` and ``--fn`` into single positional argument ``bams`` for following commands: :command:`compute-control-statistics`, :command:`compute-target-depth`, :command:`prepare-depth-of-coverage`.
* Rename ``output`` argument to ``copy-number`` for :command:`compute-copy-number` command.
* Rename ``output`` argument to ``read-depth`` for :command:`compute-read-depth` command.
* Combine optional arguments ``--gene`` and ``--region`` into single positional argument ``gene`` for :command:`compute-control-statistics` command.
* Deprecate :meth:`sdk.utils.parse_input_bams` method.
* Update :meth:`api.utils.predict_alleles` method to match ``0.31.0`` version of ``fuc`` package.
* Fix bug in :command:`filter-samples` command when ``--exclude`` argument is used for archive files with SampleTable type.
* Improve CNV caller for CYP2A6, CYP2B6, CYP2D6, CYP2E1, CYP4F2, GSTM1, SLC22A2, SULT1A1, UGT1A4, UGT2B15, and UGT2B17.
* Add a new CNV call for CYP2D6: ``PseudogeneDeletion``.
* In CYP2E1 CNV nomenclature, ``PartialDuplication`` has been renamed to ``PartialDuplicationHet`` and a new CNV call ``PartialDuplicationHom`` has been added. Furthermore, calling algorithm for CYP2E1\*S1 allele has been updated. When partial duplication is present, from now on the algorithm requires only \*7 to call \*S1 instead of both \*7 and \*4.
* Add a new CNV call for SLC22A2: ``Intron9Deletion,Exon11Deletion``.
* Add a new CNV call for UGT1A4: ``Intron1PartialDup``.
* Add new CNV calls for UGT2B15: ``PartialDeletion3`` and ``Deletion``.
* Add a new CNV call for UGT2B17: ``Deletion,PartialDeletion2``. Additionally, several CNV calls have been renamed: ``Normal`` → ``Normal,Normal``; ``DeletionHet`` → ``Normal,Deletion``; ``DeletionHom`` → ``Deletion,Deletion``; ``PartialDeletionHet`` → ``Deletion,PartialDeletion1``.

0.12.0 (2022-01-29)
-------------------

* Update :command:`run-ngs-pipeline` command to allow users to provide a custom CNV caller.
* Update :meth:`api.core.predict_phenotype` method to not raise an error when a given star allele does not exist in the allele table. From now on, the method will output a warning about it but still produce an ``Indeterminate`` call.
* Fix minor bug with ``--samples`` argument in commands :command:`plot-bam-copy-number`, :command:`plot-bam-read-depth`, :command:`plot-vcf-allele-fraction`, and :command:`plot-vcf-read-depth`.
* Update :meth:`sdk.utils.add_cn_samples` method to accept a list of samples in addition to a file.
* Add new argument ``--fontsize`` to :command:`plot-bam-read-depth` command.
* Fix minor bug in :command:`plot-bam-read-depth` command.
* Moved 1KGP reference haplotype panels and CNV callers to the ``pypgx-bundle`` `repository <https://github.com/sbslee/pypgx-bundle>`__ (only those files were moved; other files such as ``allele-table.csv`` and ``variant-table.csv`` are intact). From now on, the user needs to clone the ``pypgx-bundle`` repository with matching PyPGx version to their home directory in order for PyPGx to correctly access the moved files. This is undoubtedly annoying, but absolutely necessary for portability reasons because PyPGx has been growing exponentially in file size due to the increasing number of genes supported and their CNV complexity, to the point where it now exceeds upload size limit for PyPI (100 Mb). After removal of those files, the size of PyPGx has reduced from >100 Mb to <1 Mb.
* Add CNV caller for G6PD (mostly for sex determination since it's located on X chromosome).
* Improve CNV caller for CYP2A6, CYP2B6, CYP2D6, CYP2E1, GSTM1, SULT1A1, UGT2B15, and UGT2B17.
* Add new CNV calls for CYP2A6: ``Duplication2``, ``Duplication3``, ``Deletion2Het``, ``Deletion3Het``, ``PseudogeneDuplication``, ``Hybrid2``, ``Hybrid3``. Additionally, some CNV calls have been renamed: ``Hybrid`` → ``Hybrid1``; ``Duplication`` → ``Duplication1``; ``DeletionHet`` → ``Deletion1Het``; ``DeletionHom`` → ``Deletion1Hom``.
* Add a new CNV call for CYP2B6: ``Duplication``.
* Add new CNV calls for CYP2D6: ``Unknown1``, ``Tandem1B``, ``Multiplication``. Additionally, some CNV calls have been renamed: ``Tandem1`` → ``Tandem1A``; ``DeletionHet,Tandem1`` → ``DeletionHet,Tandem1A``; ``Duplication,Tandem1`` → ``Duplication,Tandem1A``.
* Add a new CNV call for CYP2E1: ``Duplication2``. Additionally, a CNV call have been renamed: ``Duplication`` → ``Duplication1``.
* Add new CNV calls for GSTM1: ``UpstreamDeletionHet`` and ``DeletionHet,UpstreamDeletionHet``.
* Add a new CNV call for UGT2B15: ``PartialDeletion2``. Additionally, a CNV call have been renamed: ``PartialDeletion`` → ``PartialDeletion1``.
* Add a new CNV call for UGT2B17: ``PartialDeletionHet``.

0.11.0 (2022-01-01)
-------------------

* Fix minor bug in :command:`compute-copy-number` command.
* Update :command:`plot-cn-af` command to check input files more rigorously.
* Add new method :meth:`sdk.utils.add_cn_samples`.
* Update :command:`compare-genotypes` command to output CNV comparisonw results as well.
* Update :command:`estimate-phase-beagle` command. From now on, the 'chr' prefix in contig names (e.g. 'chr1' vs. '1') will be automatically added or removed as necessary to match the reference VCF’s contig names.
* Add index files for 1KGP reference haplotype panels.
* Add new argument ``--panel`` to :command:`run-chip-pipeline` command.
* Remove 1KGP reference haplotype panels for GSTT1 and UGT2B17 because these genes only have star alleles defined with SV.
* Change 1KGP reference haplotype panels for GRCh38. Previously, PyPGx was using the panels from `Lowy-Gallego et al., 2019 <https://wellcomeopenresearch.org/articles/4-50>`__ where the authors had aligned sequence reads against the full GRCh38 reference, including ALT contigs, decoy, and EBV/IMGT/HLA sequences. This resulted in poor phasing/imputation performance for highly polymorphic PGx genes (e.g. CYP2D6) presumably because the panels were missing haplotype information for lots of SNVs/indels as sequence reads with those variants were mapped to ALT contigs; however, the panels were still the best option at the time (definitely better than lifting over GRCh37 panels). Fortunately, `Byrska-Bishop et al., 2021 <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__ from New York Genome Center has recently published a new set of GRCh38 panels which apparently has less of this problem despite still having sequence reads aligned in the presence of ALT contigs, etc. When empirically tested, these panels showed a significant increase in phasing/imputation performance. Therefore, from now on, PyPGx will use these panels for GRCh38 data.
* Update GRCh38 variant information for following alleles: CYP2D6\*35, CYP2D6\*45, CYP2D6\*46.
* Update gene region for SLC22A2 to match GRCh37 and GRCh38.
* Add CNV caller for CYP4F2 and SULT1A1.
* Improve CNV caller for CYP2A6, CYP2D6, and SLC22A2.
* Add a new CNV call for CYP2D6: ``Tandem3``.

0.10.1 (2021-12-20)
-------------------

* Fix major bug where CNV callers are not packaged properly.

0.10.0 (2021-12-19)
-------------------

* :issue:`32`: Update :command:`import-variants` command to accept phased VCF as input. It will output VcfFrame[Consolidated] if the input VCF is fully phased or otherwise VcfFrame[Imported] as usual.
* Add new property ``sdk.utils.Archive.type`` to quickly access the archive's semantic type.
* Update :meth:`sdk.utils.Archive.check_type` method to be able to test more than one semantic type at once.
* Update :meth:`api.plot.plot_vcf_allele_fraction` method to accept both VcfFrame[Imported] and VcfFrame[Consolidated].
* :issue:`32`: Update :command:`run-ngs-pipeline` command to accept phased VCF as input. In this case, the command will skip statistical haplotype phasing.
* :issue:`34`: Update commands :command:`run-ngs-pipeline` and :command:`run-chip-pipeline` to load large VCF files significantly faster by allowing random access. This also means, from now on, input VCF files must be BGZF compressed (.gz) and indexed (.tbi).
* :issue:`36`: Update phenotype data for CACNA1S, CFTR, IFNL3, RYR1 (thanks `@NTNguyen13 <https://github.com/NTNguyen13>`__).
* :pr:`39`: Add new gene F5 (thanks `@NTNguyen13 <https://github.com/NTNguyen13>`__).
* Update :command:`import-variants` command to be able to subset/exclude specified samples.
* Update :command:`import-read-depth` command to be able to subset/exclude specified samples.
* Rename ``--samples`` argument from :command:`compute-copy-number` command to ``--samples-without-sv``.
* Rename ``--samples`` argument from :command:`run-ngs-pipeline` command to ``--samples-without-sv``.
* Update :command:`run-ngs-pipeline` and :command:`run-chip-pipeline` commands to be able to subset/exclude specified samples.
* Remove ``--fn`` argument from :command:`filter-samples` command.
* Update :meth:`api.plot.plot_cn_af` method to accept both VcfFrame[Imported] and VcfFrame[Consolidated].
* Improve CNV caller for CYP2D6, GSTM1, and UGT1A4.
* Add a new CNV call for CYP2D6: ``Tandem2C``, ``DeletionHom``.
* Add a new CNV call for UGT1A4: ``Intron1DeletionB``. Additionally, a CNV call have been renamed: ``Intron1Deletion`` → ``Intron1DeletionA``.

0.9.0 (2021-12-05)
------------------

* Add 1KGP reference haplotype panels for GRCh37 for the 17 recently added genes (in ``v0.8.0``).
* Add 1KGP reference haplotype panels for GRCh38 for all target genes.
* Add GRCh37 CNV caller for UGT1A4.
* Add GRCh38 CNV callers for all ten SV genes (CYP2A6, CYP2B6, CYP2D6, CYP2E1, GSTTM1, GSTT1, SLC22A2, UGT1A4, UGT2B15, UGT2B17).
* Update commands :command:`run-ngs-pipeline` and :command:`run-chip-pipeline` to support GRCh38.
* Update the **GeT-RM WGS** tutorial to include a non-SV gene (i.e. CYP3A5) and also GRCh38.
* Make the profiles (e.g. copy number) look prettier.
* Rename :meth:`sdk.utils.Archive.check` method to :meth:`sdk.utils.Archive.check_type`.
* Add new method :meth:`sdk.utils.Archive.check_metadata`.
* Add new error ``sdk.utils.IncorrectMetadataError``.
* Update :command:`run-ngs-pipeline` command to check input files more vigorously.
* Add new method :meth:`sdk.utils.compare_metadata`.
* Add new method :meth:`api.core.get_strand`.
* Add new method :meth:`api.core.get_exon_starts`.
* Add new method :meth:`api.core.get_exon_ends`.
* :pr:`31`: Fix minor bug in commands :command:`run-ngs-pipeline` and :command:`import-read-depth` (thanks `@NTNguyen13 <https://github.com/NTNguyen13>`__).
* Fix minor bug in :meth:`api.core.predict_score` method.
* Update variant information for following alleles: CYP2D6\*27, CYP2D6\*32, CYP2D6\*131, CYP2D6\*141.

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
* Add new argument ``--fontsize`` to :command:`plot-bam-copy-number` command.
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
