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
        bam2gt2   convert BAM files to genotype files [SGE]
        gt2pt     convert a genotype file to phenotypes
        bam2vcf   convert BAM files to a VCF file
        bam2vcf2  convert BAM files to a VCF file [SGE]
        bam2gdf   convert BAM files to a GDF file
        gt2html   convert a genotype file to an HTML report
        bam2html  convert a BAM file to an HTML report [SGE]
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
        compgt    compute the concordance between two genotype files
        compvcf   calculate the concordance between two VCF files

    optional arguments:
      -h, --help  show this help message and exit
      --version   print the PyPGx version number and exit

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
      --bam_dir DIR    use all BAM files in this directory as input
      --bam_list FILE  list of input BAM files, one file per line

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
For running within Python:

>>> from pypgx.gt2pt import phenotyper
>>> phenotyper("cyp2d6", "*1", "*1")
'normal_metabolizer'
>>> phenotyper("cyp2d6", "*1", "*4")
'intermediate_metabolizer'
>>> phenotyper("cyp2d6", "*1", "*2x2")
'ultrarapid_metabolizer'
>>> phenotyper("cyp2d6", "*5", "*2x2")
'normal_metabolizer'
