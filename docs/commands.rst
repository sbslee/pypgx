Commands
********

This section describes how to use PyPGx as a command-line program.
For the public API of Python module ``pypgx``, please see the API section.

pgkb command
============

Synopsis
--------

pypgx pgkb *[options]*

Description
-----------

Extract CPIC guidelines using PharmGKB API. There are no required 
arguments.

Options
-------

-o FILE     output to FILE [stdout]
-t          extract first three guidelines for testing

report command
==============

Synopsis
--------

pypgx report *[options] gt*

Description
-----------

Create HTML report using Stargazer data. *gt* is genotype file.

Options
-------

-o FILE     output to FILE [stdout]

sdf2gdf command
===============

Synopsis
--------

pypgx sdf2gdf *[options] sdf id [id ...]*

Description
-----------

Create GDF file from SDF file. *sdf* is SDF file. *id* is sample ID.

Options
-------

-o FILE     output to FILE [stdout]

bam2sdf command
===============

Synopsis
--------

pypgx bam2sdf *[options] gb tg cg bam [bam ...]*

Description
-----------

Create SDF file from BAM file(s). *gb* is genome build (hg19, hg38). 
*tg* is target gene. *cg* is control gene or region. *bam* is BAM file.

Options
-------

-o FILE     output to FILE [stdout]

bam2gdf command
===============

Synopsis
--------

pypgx bam2gdf *[options] gb tg cg bam [bam ...]*

Description
-----------

Create GDF file from BAM file(s). *gb* is genome build (hg19, hg38). 
*tg* is target gene. *cg* is control gene or region. *bam* is BAM file.

Options
-------

-o FILE     output to FILE [stdout]

minivcf command
===============

Synopsis
--------

pypgx minivcf *[options] vcf region*

Description
-----------

Slice VCF file. *vcf* is VCF file. *region* is target region.

Options
-------

-o FILE     output to FILE [stdout]

merge command
=============

Synopsis
--------

pypgx merge *[options] vcf [vcf ...]*

Description
-----------

Merge VCF files. *vcf* is VCF file.

Options
-------

-r STR      target region
-o FILE     output to FILE [stdout]

summary command
===============

Synopsis
--------

pypgx summary *[options] tg gt*

Description
-----------

Create summary file using Stargazer data. *tg* is target gene. 
*gt* is genotype file.

Options
-------

-o FILE     output to FILE [stdout]

meta command
============

Synopsis
--------

pypgx meta *[options] tg sf [sf ...]*

Description
-----------

Create meta file from summary files. *tg* is target gene. *sf* is 
summary file.

Options
-------

-o FILE     output to FILE [stdout]

compare command
===============

Synopsis
--------

pypgx compare *[options] gt [gt ...]*

Description
-----------

Compare genotype files. *gt* is genotype file.

Options
-------

-o FILE     output to FILE [stdout]

remap command
=============

Synopsis
--------

pypgx remap *[options] conf*

Description
-----------

Remap BAM file(s) to different reference. *conf* is configuration file. 
See the API section for details.

Options
-------

There are no options.

fq2bam command
==============

Synopsis
--------

pypgx fq2bam *[options] conf*

Description
-----------

Create BAM file(s) from FASTQ file(s). *conf* is configuration file. 
See the API section for details.

Options
-------

There are no options.

sges command
============

Synopsis
--------

pypgx sges *[options] conf*

Description
-----------

Run per-sample genotyping with Stargazer. *conf* is configuration file. 
See the API section for details.

Options
-------

There are no options.

sgep command
============

Synopsis
--------

pypgx sgep *[options] conf*

Description
-----------

Run per-project genotyping with Stargazer (1). *conf* is configuration file. 
See the API section for details.

Options
-------

There are no options.

sgea command
============

Synopsis
--------

pypgx sgea *[options] conf*

Description
-----------

Run per-project genotyping with Stargazer (2). *conf* is configuration file. 
See the API section for details.

Options
-------

There are no options.

cpa command
===========

Synopsis
--------

pypgx cpa *[options] rdata*

Description
-----------

Run change point analysis for copy number. *rdata* is Rdata file.

Options
-------

-o FILE     output to FILE [stdout]

plotcov command
===============

Synopsis
--------

pypgx plotcov *[options] sdf out*

Description
-----------

Plot coverage data to PDF file. *sdf* is SDF file. *out* is PDF file.

Options
-------

There are no options.

check command
=============

Synopsis
--------

pypgx check *[options] star snp*

Description
-----------

Check table files for Stargazer. *star* is star allele table file. 
*snp* is SNP table file.

Options
-------

There are no options.

liftover command
================

Synopsis
--------

pypgx liftover *[options] star snp tg*

Description
-----------

Convert variants in SNP table from hg19 to hg38. *star* is star allele
table file. *snp* is SNP table file. *tg* is target gene.

Options
-------

-o FILE     output to FILE [stdout]

peek command
============

Synopsis
--------

pypgx peek *[options] vcf*

Description
-----------

Find all possible star alleles from VCF file. *vcf* is VCF file.

Options
-------

-o FILE     output to FILE [stdout]

snp command
===========

Synopsis
--------

pypgx snp *[options] vcf pair [pair ...]*

Description
-----------

View variant data for sample/star allele pairs. *vcf* is VCF file. 
*pair* is sample/star allele pair.

Options
-------

-o FILE     output to FILE [stdout]

bam2vcf command
===============

Synopsis
--------

pypgx bam2vcf *[options] gb tg fa bam [bam ...]*

Description
-----------

Create a single- or multi-sample VCF file from one or more BAM files.
*gb* is the genome build ('hg19' or 'hg38'). *tg* is the target gene 
(e.g. 'cyp2d6') or region (e.g. 'chr22:42512500-42551883'). 
*fa* is the reference FASTA file. *bam* is the input BAM file(s). 

.. note::

    This is essentially a wrapper for the ``mpileup`` and ``call`` 
    commands from the BCFtools program; therefore, you must have 
    BCFtools installed before running this command.

Options
-------

-o FILE     output to FILE [stdout]