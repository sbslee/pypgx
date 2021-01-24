README
******

.. image:: https://badge.fury.io/py/pypgx.svg
    :target: https://badge.fury.io/py/pypgx
.. image:: https://readthedocs.org/projects/pypgx/badge/?version=latest
    :target: https://pypgx.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

Table of Contents
=================

* `Introduction`_
* `Installation`_
* `pypgx CLI`_
* `pypgx API`_

Introduction
============

pypgx is a Python package for pharmacogenomics research, which can be used as a standalone program and as a Python module. Documentation is available at `Read the Docs <https://pypgx.readthedocs.io/en/latest/>`_.

Installation
============

You can easily install pypgx and all of its dependencies with the Anaconda distribution.

.. code-block:: console

   $ conda create -n pypgx -c sbslee pypgx

Before using pypgx, make sure to activate the conda environment where pypgx is installed.

.. code-block:: console

  $ conda activate pypgx

pypgx CLI
=========

The `pypgx CLI page <https://pypgx.readthedocs.io/en/latest/cli.html>`_ describes command-line interface (CLI) for the pypgx package.

You can display help message for pypgx CLI by entering:

.. code-block:: console

    $ pypgx -h
    usage: pypgx [-v] [-h] COMMAND ...

    positional arguments:
      COMMAND               Name of the command.
        compare-stargazer-calls
                            Compute the concordance between two 'genotype-
                            calls.tsv' files created by Stargazer.
        calculate-read-depth
                            Create a GDF (GATK DepthOfCoverage Format) file for
                            Stargazer from BAM files by computing read depth.
        call-variants-gatk-sge
                            Create a VCF (Variant Call Format) file for Stargazer
                            from BAM files by calling SNVs and indels.

    optional arguments:
      -v, --version         Show the version and exit.
      -h, --help            Show this help message and exit.


You can display command-specific help message by entering (e.g. ``calculate-read-depth``):

.. code-block:: console

    $ pypgx calculate-read-depth -h
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
      -h, --help            Show this help message and exit.

For running in command line:

.. code-block:: console

    $ pypgx calculate-read-depth -t cyp2d6 -c vdr -i bam-list.txt -o read-depth.gdf

The output GDF file will look something like:

.. parsed-literal::

    Locus	Total_Depth	Average_Depth_sample	Depth_for_Steven	Depth_for_John
    ...
    chr22:42539471	190	95	53	137
    chr22:42539472	192	96	54	138
    chr22:42539473	190	95	53	137
    ...

pypgx API
=========

The `pypgx API page <https://pypgx.readthedocs.io/en/latest/api.html>`_ describes application programming interface (API) for the pypgx package.

For running within Python (e.g. ``phenotyper``):

.. code:: ipython3

    from pypgx import phenotyper
    print(phenotyper("cyp2d6", "*1", "*1"))
    print(phenotyper("cyp2d6", "*1", "*4"))
    print(phenotyper("cyp2d6", "*1", "*2x2"))  # *2x2 is gene duplication.
    print(phenotyper("cyp2d6", "*4", "*5"))    # *5 is gene deletion.

To give:

.. parsed-literal::

    normal_metabolizer
    intermediate_metabolizer
    ultrarapid_metabolizer
    poor_metabolizer
