PyPGx
*****

.. image:: https://badge.fury.io/py/pypgx.svg
    :target: https://badge.fury.io/py/pypgx

PyPGx is a Python package for pharmacogenomics (PGx) research.

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

    optional arguments:
      -h, --help  show this help message and exit

For getting tool-specific help (e.g. ``pgkb`` tool)::

    $ pypgx pgkb -h
    usage: pypgx pgkb [-h] [-o FILE]

    optional arguments:
      -h, --help  show this help message and exit
      -o FILE     output to FILE [stdout]

For running in command line (e.g. ``pgkb`` tool)::

    $ pypgx pgkb -o guidelines.txt

Running within Python
=====================
For running within Python (e.g. ``pgkb`` tool):

>>> from pypgx.pgkb import pgkb
>>> guidelines = pgkb("pgkb")
>>> print(guidelines)
