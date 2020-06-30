Commands
********

This section describes how to use PyPGx as a command-line program.
For the public API of Python module ``pypgx``, please see the API section.

pgkb
====

Synopsis
--------

pypgx pgkb [options]

Description
-----------

Extract CPIC guidelines using PharmGKB API. There are no required 
arguments.

Options
-------

-o FILE     output to FILE [stdout]
-t          will only look first three guidelines

report
======

Synopsis
--------

pypgx report [options] gt

Description
-----------

Create HTML report using data from Stargazer. *gt* is genotype file 
from Stargazer.

Options
-------

-o FILE     output to FILE [stdout]

sdf2gdf
=======

Synopsis
--------

pypgx sdf2gdf [options] sdf id [id ...]

Description
-----------

Create GDF file from SDF file. *sdf* is SDF file. *id* is sample ID.

Options
-------

-o FILE     output to FILE [stdout]

bam2sdf
=======

Synopsis
--------

pypgx bam2sdf [options] tg cg bam [bam ...]

Description
-----------

Create SDF file from BAM file(s). *tg* is target gene. *cg* is control 
gene. *bam* is BAM file.

Options
-------

-o FILE     output to FILE [stdout]

bam2gdf
=======

Synopsis
--------

pypgx bam2gdf [options] tg cg bam [bam ...]

Description
-----------

Create GDF file from BAM file(s). *tg* is target gene. *cg* is control 
gene. *bam* is BAM file.

Options
-------

-o FILE     output to FILE [stdout]

minivcf
=======

Synopsis
--------

pypgx minivcf [options] vcf region

Description
-----------

Slice VCF file. *vcf* is VCF file. *region* is target region.

Options
-------

-o FILE     output to FILE [stdout]

merge
========

Synopsis
--------

pypgx merge [options] vcf [vcf ...]

Description
-----------

Merge VCF files. *vcf* is VCF file.

Options
-------

-r STR      target region
-o FILE     output to FILE [stdout]

summary
=======

Synopsis
--------

pypgx summary [options] gt tg

Description
-----------

Create summary file using data from Stargazer. *gt* is genotype file 
from Stargazer. *tg* is target gene.

Options
-------

-o FILE     output to FILE [stdout]

meta
====

Synopsis
--------

pypgx meta [options] tg sf [sf ...]

Description
-----------

Create meta file from summary files. *tg* is target gene. *sf* is 
summary file.

Options
-------

-o FILE     output to FILE [stdout]