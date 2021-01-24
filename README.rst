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
* `Running in Command Line`_
* `Running within Python`_

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

Running in Command Line
=======================

For getting help, enter ``pypgx -h``:

.. code-block:: console

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



For running in command line::

    $ pypgx bam2gdf hg19 cyp2d6 vdr out.gdf in1.bam in2.bam

The output GDF file will look like::

    Locus	Total_Depth	Average_Depth_sample	Depth_for_S1	Depth_for_S2
    ...
    chr22:42539471	190	95	53	137
    chr22:42539472	192	96	54	138
    chr22:42539473	190	95	53	137
    ...

Running within Python
=====================

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
