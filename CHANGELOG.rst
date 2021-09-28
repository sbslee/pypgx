Changelog
*********

0.5.0 (in development)
----------------------

* Update ``api.utils.create_read_depth_tsv`` method to automatically detect ``chr`` string in input BAM.
* Add ``sdk.utils.parse_input_bams`` method.
* Add the 1000 Genomes Project reference haplotype panel for GRCh37. When estimating haplotype phase of observed variants, users are no longer needed to download and specify a panel. GRCh38 support will follow in a future release.
* Rename command :command:`create-read-depth-tsv` to :command:`prepare-depth-of-coverage`.
* Add ``bed`` argument to :command:`prepare-depth-of-coverage` command.
* Update :command:`prepare-depth-of-coverage` command to output archive file instead of TSV file.
* Update :command:`import-read-depth` command to accept archive file as input instead of TSV file.
* Add ``fitted`` argument to ``api.plot.plot_bam_copy_number`` method.
* From now on, missing copy number will be imputed with forward filling instead of column median.
* Update :command:`predict-cnv` command to support a user-defined CNV caller.

0.4.1 (2021-09-21)
------------------

* Initial release.
