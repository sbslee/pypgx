Changelog
*********

0.5.0 (in development)
----------------------

* Update ``api.utils.create_read_depth_tsv`` method to automatically detect ``chr`` string in input BAM.
* Add ``sdk.utils.parse_input_bams`` method.
* Add the 1000 Genomes Project reference haplotype panel for GRCh37. When estimating haplotype phase of observed variants, users are no longer needed to download and specify a panel. GRCh38 support will follow in a future release.
* Rename command :command:`create-read-depth-tsv` to :command:`prepare-depth-of-coverage`.

0.4.1 (2021-09-21)
------------------

* Initial release.
