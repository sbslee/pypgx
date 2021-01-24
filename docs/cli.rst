CLI
***

This section describes command-line interface (CLI) for the pypgx package.

compare-stargazer-calls
=======================

.. code-block:: console

    usage: pypgx compare-stargazer-calls [-h] -r PATH -t PATH -o PATH

    Compute the concordance between two 'genotype-calls.tsv' files created by
    Stargazer.

    optional arguments:
      -h, --help            show this help message and exit
      -r PATH, --ref-file PATH
                            Path to the reference or truth 'genotype-calls.tsv'
                            file created by Stargazer. [required]
      -t PATH, --test-file PATH
                            Path to the test 'genotype-calls.tsv' file created by
                            Stargazer. [required]
      -o PATH, --output-file PATH
                            Path to the output file. [required]

calculate-read-depth
====================

.. code-block:: console

    usage: pypgx calculate-read-depth -t TEXT -c TEXT [-i PATH] -o PATH [-a TEXT]
                                      [-h]

    Create a GDF (GATK DepthOfCoverage Format) file for Stargazer from BAM files
    by computing read depth.

    optional arguments:
      -t TEXT, --target-gene TEXT
                            Name of the target gene. Choices: {'abcb1', 'cacna1s',
                            'cftr', 'cyp1a1', 'cyp1a2', 'cyp1b1', 'cyp2a6',
                            'cyp2a13', 'cyp2b6', 'cyp2c8', 'cyp2c9', 'cyp2c19',
                            'cyp2d6', 'cyp2e1', 'cyp2f1', 'cyp2j2', 'cyp2r1',
                            'cyp2s1', 'cyp2w1', 'cyp3a4', 'cyp3a5', 'cyp3a7',
                            'cyp3a43', 'cyp4a11', 'cyp4a22', 'cyp4b1', 'cyp4f2',
                            'cyp17a1', 'cyp19a1', 'cyp26a1', 'dpyd', 'g6pd',
                            'gstm1', 'gstp1', 'gstt1', 'ifnl3', 'nat1', 'nat2',
                            'nudt15', 'por', 'ptgis', 'ryr1', 'slc15a2',
                            'slc22a2', 'slco1b1', 'slco1b3', 'slco2b1', 'sult1a1',
                            'tbxas1', 'tpmt', 'ugt1a1', 'ugt1a4', 'ugt2b7',
                            'ugt2b15', 'ugt2b17', 'vkorc1', 'xpc'}. [required]
      -c TEXT, --control-gene TEXT
                            Name of a preselected control gene. Used for
                            intrasample normalization during copy number analysis
                            by Stargazer. Choices: {'egfr', 'ryr1', 'vdr'}.
                            Alternatively, you can provide a custom genomic region
                            with the 'chr:start-end' format (e.g.
                            chr12:48232319-48301814). [required]
      -i PATH, --bam-path PATH
                            Read BAM files from PATH, one file path per line.
                            [required]
      -o PATH, --output-file PATH
                            Path to the output file. [required]
      -a TEXT, --genome-build TEXT
                            Build of the reference genome assembly. Choices:
                            {'hg19', 'hg38'}. [default: 'hg19']
      -h, --help                Show this help message and exit.
