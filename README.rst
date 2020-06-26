PyPGx
*****

.. image:: https://badge.fury.io/py/pypgx.svg
    :target: https://badge.fury.io/py/pypgx

PyPGx is a Python package for pharmacogenomics (PGx) research. Documentation is available at `Read the Docs <https://pypgx.readthedocs.io/en/latest/>`_.

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
>>> gdf = bam2gdf("bam2gdf", "cyp2d6", "vdr", "in1.bam", "in2.bam", "in3.bam")
>>> for line in gdf.split("\n"):
>>>     print(line)

To give::

    Locus	Total_Depth	Average_Depth_sample	Depth_for_S1	Depth_for_S2	Depth_for_S3
    ...
    chr22:42526563	858	286	237	432	189
    chr22:42526564	860	286.67	239	433	188
    chr22:42526565	857	285.67	239	433	185
    ...
