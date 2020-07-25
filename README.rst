PyPGx
*****

.. image:: https://badge.fury.io/py/pypgx.svg
    :target: https://badge.fury.io/py/pypgx
.. image:: https://readthedocs.org/projects/pypgx/badge/?version=latest
    :target: https://pypgx.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

PyPGx is a Python package for pharmacogenomics (PGx) research, which can be 
used as a standalone program and as a Python module. Documentation is 
available at `Read the Docs <https://pypgx.readthedocs.io/en/latest/>`_.

Stargazer
=========

For genotype analyses PyPGx relies on Stargazer, a bioinformatics tool for 
calling star alleles (haplotypes) in PGx genes using data from 
next-generation sequencing (NGS) or single nucleotide polymorphism (SNP) 
array. For more information, please visit Stargazer's 
`official webpage <https://stargazer.gs.washington.edu/stargazerweb>`_ and 
`Github repository <https://github.com/sbslee/stargazer>`_.

Installation
============

The easiest way to install PyPGx is to use ``pip``::

    $ pip install pypgx

Running in Command Line
=======================

For getting help::

    $ pypgx -h
    usage: pypgx [-h] [--version] tool ...

    positional arguments:
      tool        name of the tool
        bam2gt    convert BAM files to a genotype file
        gt2pt     convert a genotype file to phenotypes
        bam2vcf   convert BAM files to a VCF file
        bam2gdf   convert BAM files to a GDF file
        gt2html   convert a genotype file to an HTML report
        bam2html  convert a BAM file to an HTML report [SGE]
        sgep      convert BAM files to a genotype file [SGE]
        xgep      convert BAM files to genotype files [SGE]
        fq2bam    convert FASTQ files to BAM files [SGE]
        bam2bam   realign BAM files to another reference genome [SGE]
        bam2sdf   convert BAM files to a SDF file
        sdf2gdf   convert a SDF file to a GDF file
        pgkb      extract CPIC guidelines using PharmGKB API
        minivcf   slice VCF file
        merge     merge VCF files
        summary   create summary file using Stargazer data
        meta      create meta file from summary files
        compare   compare genotype files
        cpa       run change point analysis for copy number
        plotcov   plot coverage data to PDF file
        check     check table files for Stargazer
        liftover  convert variants in SNP table from hg19 to hg38
        peek      find all possible star alleles from VCF file
        snp       view variant data for sample/star allele pairs

    optional arguments:
      -h, --help  show this help message and exit
      --version   print the PyPGx version number and exit

For getting tool-specific help::

    $ pypgx bam2gdf -h
    usage: pypgx bam2gdf [-h] [--bd DIR] [--bl FILE] [-o FILE]
                         gb tg cg [bam [bam ...]]

    positional arguments:
      gb          genome build
      tg          target gene
      cg          control gene or region
      bam         BAM file

    optional arguments:
      -h, --help  show this help message and exit
      --bd DIR    directory containing BAM files
      --bl FILE   list of BAM files, one file per line
      -o FILE     output to FILE [stdout]

For running in command line::

    $ pypgx bam2gdf hg19 cyp2d6 vdr in1.bam in2.bam -o out.gdf

Running within Python
=====================
For running within Python:

>>> from pypgx.bam2gdf import bam2gdf
>>> bam = ["in1.bam", "in2.bam"]
>>> gdf = bam2gdf("hg19", "cyp2d6", "vdr", bam)
>>> for line in gdf.split("\n"):
>>>     print(line)

To give::

    Locus	Total_Depth	Average_Depth_sample	Depth_for_S1	Depth_for_S2
    ...
    chr22:42539471	190	95	53	137
    chr22:42539472	192	96	54	138
    chr22:42539473	190	95	53	137
    ...
