pypgx
*****

.. image:: https://badge.fury.io/py/pypgx.svg
    :target: https://badge.fury.io/py/pypgx
.. image:: https://readthedocs.org/projects/pypgx/badge/?version=latest
    :target: https://pypgx.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

Table of Contents
=================

* `Introduction`_
* `Installation`_
* `Stargazer`_
* `Sun Grid Engine (SGE)`_
* `SNP Callers`_
* `pypgx CLI`_
* `pypgx API`_

Introduction
============

pypgx is a Python package for pharmacogenomics research, which can be used as a standalone program and as a Python module. Documentation is available at `Read the Docs <https://pypgx.readthedocs.io/en/latest/>`_.

Installation
============

You can easily install pypgx and all of its dependencies with the Anaconda distribution.

.. code-block:: console

   conda create -n pypgx -c sbslee pypgx
   conda activate pypgx

Stargazer
=========

For genotype analyses pypgx relies on Stargazer, a bioinformatics tool for
calling star alleles (haplotypes) in PGx genes using data from
next-generation sequencing (NGS) or single nucleotide polymorphism (SNP)
array. Therefore, Stargazer must be pre-installed in order to run pypgx
commands such as ``bam2gt``. For more information on Stargazer, please visit
their `official webpage <https://stargazer.gs.washington.edu/stargazerweb>`_
and `Github repository <https://github.com/sbslee/stargazer>`_.

Sun Grid Engine (SGE)
=====================

Many pypgx commands such as ``bam2gt2`` rely on the Sun Grid Engine (SGE)
cluster to distribute their tasks across multiple machines for speed. These
commands are indicated by ``[SGE]`` and will generate a shell script, which
can be run like this::

    $ sh example-qsub.sh

SNP Callers
===========

One major input for the Stargzer program is a Variant Call Format (VCF) file,
which is a standard file format for storing SNP calls. Currently, pypgx
relies on two SNP callers to make VCF files: Genome Analysis Toolkit (GATK)
and BCFtools. When running pypgx commands like ``bam2vcf``, you can pick
which SNP calling algorithm to use; it is assumed that you already installed
the corresponding SNP caller.

Generally speaking, GATK is considered more accurate but much slower
than BCFtools. For instance, without the use of the SGE cluster, SNP calling
for 70 WGS samples for the CYP2D6 gene takes 19 min to complete with GATK,
but only 2 min with BCFtools. Therefore, if you have many samples and you do
not have access to SGE for running parallel jobs, BCFtools may be a better
choice. Of course, if you have SGE in your sever, then GATK is strongly
recommended.

For more information on the SNP callers, please visit the
`GATK website <https://gatk.broadinstitute.org/hc/en-us>`_ and
the `BCFtools website <http://samtools.github.io/bcftools/bcftools.html>`_.

pypgx CLI
=========

You can display help message for pypgx CLI by entering::

    pypgx -h

To give::

    usage: pypgx [-v] [-h] COMMAND ...

    positional arguments:
      COMMAND               Name of the command.
        compare-stargazer-calls
                            Compute the concordance between two 'genotype-
                            calls.tsv' files created by Stargazer.
        calculate-read-depth
                            Create a GDF (GATK DepthOfCoverage Format) file for
                            Stargazer from BAM files by computing read depth.

    optional arguments:
      -v, --version         Show the version and exit.
      -h, --help            Show this help message and exit.

For getting command-specific help (e.g. ``calculate-read-depth``), enter::

    pypgx compare-stargazer-calls -h

To give::

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

For running in command line::

    pypgx calculate-read-depth \
    -t cyp2d6 \
    -c vdr \
    -i bam-list.txt \
    -o read-depth.gdf

The output GDF file will look like::

    Locus	Total_Depth	Average_Depth_sample	Depth_for_Steven	Depth_for_John
    ...
    chr22:42539471	190	95	53	137
    chr22:42539472	192	96	54	138
    chr22:42539473	190	95	53	137
    ...

pypgx API
=========

.. code:: ipython3

    results = pypgx.calculate_read_depth("cyp2d6", "vdr", "bam-list.txt")

.. code:: ipython3

    results.df

    .. raw:: html

        <div>
        <style scoped>
            .dataframe tbody tr th:only-of-type {
                vertical-align: middle;
            }

            .dataframe tbody tr th {
                vertical-align: top;
            }

            .dataframe thead th {
                text-align: right;
            }
        </style>
        <table border="1" class="dataframe">
          <thead>
            <tr style="text-align: right;">
              <th></th>
              <th>Locus</th>
              <th>Total_Depth</th>
              <th>Average_Depth_sample</th>
              <th>Depth_for_74608</th>
              <th>Depth_for_74608</th>
            </tr>
          </thead>
          <tbody>
            <tr>
              <th>0</th>
              <td>12:48232319</td>
              <td>0</td>
              <td>0.0</td>
              <td>0</td>
              <td>0</td>
            </tr>
            <tr>
              <th>1</th>
              <td>12:48232320</td>
              <td>0</td>
              <td>0.0</td>
              <td>0</td>
              <td>0</td>
            </tr>
            <tr>
              <th>2</th>
              <td>12:48232321</td>
              <td>0</td>
              <td>0.0</td>
              <td>0</td>
              <td>0</td>
            </tr>
            <tr>
              <th>3</th>
              <td>12:48232322</td>
              <td>0</td>
              <td>0.0</td>
              <td>0</td>
              <td>0</td>
            </tr>
            <tr>
              <th>4</th>
              <td>12:48232323</td>
              <td>0</td>
              <td>0.0</td>
              <td>0</td>
              <td>0</td>
            </tr>
            <tr>
              <th>...</th>
              <td>...</td>
              <td>...</td>
              <td>...</td>
              <td>...</td>
              <td>...</td>
            </tr>
            <tr>
              <th>108875</th>
              <td>22:42551879</td>
              <td>0</td>
              <td>0.0</td>
              <td>0</td>
              <td>0</td>
            </tr>
            <tr>
              <th>108876</th>
              <td>22:42551880</td>
              <td>0</td>
              <td>0.0</td>
              <td>0</td>
              <td>0</td>
            </tr>
            <tr>
              <th>108877</th>
              <td>22:42551881</td>
              <td>0</td>
              <td>0.0</td>
              <td>0</td>
              <td>0</td>
            </tr>
            <tr>
              <th>108878</th>
              <td>22:42551882</td>
              <td>0</td>
              <td>0.0</td>
              <td>0</td>
              <td>0</td>
            </tr>
            <tr>
              <th>108879</th>
              <td>22:42551883</td>
              <td>0</td>
              <td>0.0</td>
              <td>0</td>
              <td>0</td>
            </tr>
          </tbody>
        </table>
        <p>108880 rows × 5 columns</p>
        </div>

For running within Python::

    from pypgx.phenotyper import phenotyper
    phenotyper("cyp2d6", "*1", "*1")
    phenotyper("cyp2d6", "*1", "*4")
    phenotyper("cyp2d6", "*1", "*2x2")  # *2x2 is gene duplication.
    phenotyper("cyp2d6", "*4", "*5")    # *5 is gene deletion.

To give::

    'normal_metabolizer'
    'intermediate_metabolizer'
    'ultrarapid_metabolizer'
    'poor_metabolizer'
