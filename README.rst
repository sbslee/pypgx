PyPGx
*****

.. image:: https://badge.fury.io/py/pypgx.svg
    :target: https://badge.fury.io/py/pypgx
.. image:: https://readthedocs.org/projects/pypgx/badge/?version=latest
    :target: https://pypgx.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

PyPGx is a Python package for pharmacogenomics (PGx) research, which can be 
used as a standalone program and as a Python module.

Documentation is available at `Read the Docs <https://pypgx.readthedocs.io/en/latest/>`_.

Installation
============

The easiest way to install PyPGx is to use ``pip``::

    $ pip install pypgx

Running in Command Line
=======================

For getting help::

    $ pypgx -h
    usage: pypgx [-h] tool ...

    positional arguments:
    tool        name of tool
      pgkb      extract CPIC guidelines using PharmGKB API
      report    create HTML report using data from Stargazer
      sdf2gdf   create GDF file from SDF file
      bam2sdf   create SDF file from BAM file(s)
      bam2gdf   create GDF file from BAM file(s)
      minivcf   slice VCF file
      merge     merge VCF files
      summary   create summary using data from Stargazer

    optional arguments:
      -h, --help  show this help message and exit

For getting tool-specific help (e.g. ``bam2gdf`` tool)::

    $ pypgx bam2gdf -h
    usage: pypgx bam2gdf [-h] [-o FILE] tg cg bam [bam ...]

    positional arguments:
      tg          target gene
      cg          control gene
      bam         BAM file

    optional arguments:
      -h, --help  show this help message and exit
      -o FILE     output to FILE [stdout]

For running in command line (e.g. ``bam2gdf`` tool)::

    $ pypgx bam2gdf -o out.gdf cyp2d6 vdr in1.bam in2.bam in3.bam

Running within Python
=====================
For running within Python (e.g. ``bam2gdf`` tool):

>>> from pypgx.bam2gdf import bam2gdf
>>> bams = ["in1.bam", "in2.bam"]
>>> gdf = bam2gdf("cyp2d6", "vdr", bams)
>>> for line in gdf.split("\n"):
>>>     print(line)

To give::

    Locus	Total_Depth	Average_Depth_sample	Depth_for_S1	Depth_for_S2
    ...
    chr22:42539471	190	95	53	137
    chr22:42539472	192	96	54	138
    chr22:42539473	190	95	53	137
    ...
