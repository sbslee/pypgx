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

For genotype analysis PyPGx relies on Stargazer, a bioinformatics tool for 
calling star alleles (haplotypes) in PGx genes using data from 
next-generation sequencing (NGS) or single nucleotide polymorphism (SNP) 
array. Stargazer is free for academic use and can be downloaded from 
`here <https://stargazer.gs.washington.edu/stargazerweb/>`_.

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
        pgkb      extract CPIC guidelines using PharmGKB API
        report    create HTML report using Stargazer data
        sdf2gdf   create GDF file from SDF file
        bam2sdf   create SDF file from BAM file(s)
        bam2gdf   create GDF file from BAM file(s)
        minivcf   slice VCF file
        merge     merge VCF files
        summary   create summary file using Stargazer data
        meta      create meta file from summary files
        compare   compare genotype files
        remap     remap BAM file(s) to different reference
        fq2bam    create BAM file(s) from FASTQ file(s)
        sges      run per-sample genotyping with Stargazer
        sgep      run per-project genotyping with Stargazer (1)
        sgea      run per-project genotyping with Stargazer (2)
        cpa       run change point analysis for copy number
        plotcov   plot coverage data to PDF file
        check     check table files for Stargazer
        liftover  convert variants in SNP table from hg19 to hg38
        peek      find all possible star alleles from VCF file
        snp       view variant data for sample/star allele pairs
        bam2vcf   create VCF file from BAM file(s)
        genotype  call star alleles from BAM file(s)
        gt2pt     call phenotypes from star alleles

    optional arguments:
      -h, --help  show this help message and exit
      --version   print the PyPGx version number and exit

For getting tool-specific help (e.g. ``bam2gdf`` tool)::

    $ pypgx bam2gdf -h
    usage: pypgx bam2gdf [-h] [--bd DIR] [--bl FILE] [-o FILE] gb tg cg [bam]

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

For running in command line (e.g. ``bam2gdf`` tool)::

    $ pypgx bam2gdf hg19 cyp2d6 vdr in1.bam in2.bam -o out.gdf

Running within Python
=====================
For running within Python (e.g. ``bam2gdf`` tool):

>>> from pypgx.bam2gdf import bam2gdf
>>> bams = ["in1.bam", "in2.bam"]
>>> gdf = bam2gdf("hg19", "cyp2d6", "vdr", bams)
>>> for line in gdf.split("\n"):
>>>     print(line)

To give::

    Locus	Total_Depth	Average_Depth_sample	Depth_for_S1	Depth_for_S2
    ...
    chr22:42539471	190	95	53	137
    chr22:42539472	192	96	54	138
    chr22:42539473	190	95	53	137
    ...
