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

Create HTML report using Stargazer data. *gt* is genotype file.

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

pypgx summary [options] tg gt

Description
-----------

Create summary file using Stargazer data. *tg* is target gene. 
*gt* is genotype file.

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

compare
=======

Synopsis
--------

pypgx compare [options] gt [gt ...]

Description
-----------

Compare genotype files. *gt* is genotype file.

Options
-------

-o FILE     output to FILE [stdout]

remap
=====

Synopsis
--------

pypgx remap [options] conf

Description
-----------

Remap BAM file(s) to different reference. *conf* is configuration file 
(see the API section for details).

Options
-------

There are no options.

fq2bam
======

Synopsis
--------

pypgx fq2bam [options] conf

Description
-----------

Create BAM file(s) from FASTQ file(s). *conf* is configuration file 
(see the API section for details).

Options
-------

There are no options.

sgea
====

Synopsis
--------

pypgx sgea [options] conf

Description
-----------

Run per-project genotyping pipeline (Anaconda-based). *conf* is 
configuration file (see the API section for details).

Options
-------

There are no options.