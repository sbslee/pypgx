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

For getting help::

    $ pypgx -h
    usage: pypgx [-h] [-v] tool ...

    positional arguments:
      tool           name of the tool
        bam2gt       convert BAM files to a genotype file
        bam2gt2      convert BAM files to genotype files [SGE]
        bam2vcf      convert BAM files to a VCF file
        bam2vcf2     convert BAM files to a VCF file [SGE]
        bam2gdf      convert BAM files to a GDF file
        gt2html      convert a genotype file to an HTML report
        bam2html     convert a BAM file to an HTML report [SGE]
        fq2bam       convert FASTQ files to BAM files [SGE]
        bam2bam      realign BAM files to another reference genome [SGE]
        sdf2gdf      convert a SDF file to a GDF file
        pgkb         extract CPIC guidelines using PharmGKB API
        minivcf      slice VCF file
        summary      create summary file using Stargazer data
        meta         create meta file from summary files
        compare      compare genotype files
        peek         find all possible star alleles from VCF file
        viewsnp      view SNP data for pairs of sample/star allele
        compvcf      compute the concordance between two VCF files
        unicov       compute the uniformity of sequencing coverage

    optional arguments:
      -h, --help     show this help message and exit
      -v, --version  print the pypgx version number and exit

For getting tool-specific help::

    $ pypgx bam2gdf -h
    usage: pypgx bam2gdf [-h] [--bam_dir DIR] [--bam_list FILE]
                         genome_build target_gene control_gene output_file
                         [bam_file [bam_file ...]]

    positional arguments:
      genome_build     genome build ('hg19' or 'hg38')
      target_gene      name of target gene (e.g. 'cyp2d6')
      control_gene     name or region of control gene (e.g. ‘vdr’,
                       ‘chr12:48232319-48301814’)
      output_file      write output to this file
      bam_file         input BAM files

    optional arguments:
      -h, --help       show this help message and exit
      --bam_dir DIR    treat any BAM files in DIR as input
      --bam_list FILE  read BAM files from FILE, one file path per line

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
