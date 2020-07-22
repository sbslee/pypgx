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

Extract CPIC guidelines using PharmGKB API.

There are no required arguments.

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

Create HTML report using Stargazer data.

*gt* is the genotype file from Stargazer.

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

Create GDF file from SDF file.

*sdf* is SDF file. *id* is sample ID.

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

Create SDF file from BAM file(s).

*gb* is the genome build ('hg19' or 'hg38'). *tg* is the target 
gene (e.g. 'cyp2d6'). *cg* is the control gene (e.g. 'vdr') or 
region (e.g. 'chr12:48232319-48301814'). *bam* is the BAM file.

Options
-------

-o FILE     output to FILE [stdout]

bam2gdf command
===============

Synopsis
--------

pypgx bam2gdf *[options] gb tg cg [bam [bam ...]]*

Description
-----------

Create GDF file from BAM file(s).

*gb* is the genome build ('hg19' or 'hg38'). *tg* is the target 
gene (e.g. 'cyp2d6'). *cg* is the control gene (e.g. 'vdr') or 
region (e.g. 'chr12:48232319-48301814'). *bam* is the input BAM file(s).

If you have many input BAM files, you may want to use the ``--bd`` or 
``--bl`` argument instead of manually listing individual files for *bam*.

Options
-------

--bd DIR    directory containing BAM files
--bl FILE   list of BAM files, one file per line
-o FILE     output to FILE [stdout]

minivcf command
===============

Synopsis
--------

pypgx minivcf *[options] vcf region*

Description
-----------

Slice VCF file.

*vcf* is VCF file. *region* is target region.

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

Merge VCF files.

*vcf* is VCF file.

Options
-------

-r STR      target region
-o FILE     output to FILE [stdout]

summary command
===============

Synopsis
--------

pypgx summary *[options] gt*

Description
-----------

Create summary file using Stargazer data.

*gt* is the genotype file from Stargazer.

Options
-------

-o FILE     output to FILE [stdout]

meta command
============

Synopsis
--------

pypgx meta *[options] sf [sf ...]*

Description
-----------

Create meta file from summary files.

*sf* is the summary file from the ``summary`` command.

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

Compare genotype files.

*gt* is the genotype file from Stargazer.

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

Remap BAM file(s) to different reference.

*conf* is the configuration file. See the API section for details.

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

Create BAM file(s) from FASTQ file(s).

*conf* is the configuration file. See the API section for details.

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

Run per-sample genotyping for multiple genes with SGE.

This command runs the per-sample genotyping pipeline by submitting 
jobs to the Sun Grid Engine (SGE) cluster. After genotype analysis by 
Stargazer, it will generate a HTML report using the ``report`` tool.

*conf* is the configuration file. See the API section for details.

.. note::

    SGE and Stargazer must be pre-installed.

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

Run per-project genotyping for single gene with SGE (1).

This command runs the per-project genotyping pipeline by submitting 
jobs to the Sun Grid Engine (SGE) cluster.

*conf* is the configuration file. See the API section for details.

.. note::

    SGE and Stargazer must be pre-installed.

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

Run per-project genotyping for single gene with SGE (2).

This command runs the per-project genotyping pipeline by submitting 
jobs to the Sun Grid Engine (SGE) cluster. The main difference between
``sgea`` and ``sgep`` is that the former uses Genome Analysis Tool 
Kit (GATK) v4 while the latter uses GATK v3.

*conf* is the configuration file. See the API section for details.

.. note::

    SGE, Stargazer, and GATK must be pre-installed.

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

Plot coverage data to PDF file.

*sdf* is SDF file. *out* is PDF file.

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

Check table files for Stargazer.

*star* is star allele table file. *snp* is SNP table file.

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

Convert variants in SNP table from hg19 to hg38.

*star* is star allele table file. *snp* is SNP table file. 
*tg* is target gene.

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

Find all possible star alleles from VCF file.

*vcf* is VCF file.

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

View variant data for sample/star allele pairs.

*vcf* is VCF file. *pair* is sample/star allele pair.

Options
-------

-o FILE     output to FILE [stdout]

bam2vcf command
===============

Synopsis
--------

pypgx bam2vcf *[options] gb tg fa [bam [bam ...]]*

Description
-----------

Create a VCF file from BAM file(s).

This command outputs a single- or multi-sample VCF file from one or 
more input BAM files. The output VCF file will only contain variants
within the target gene or region. This is essentially a wrapper with
certain parameters for various commands from the BCFtools program 
(e.g. ``mpileup`` and ``call``). This means the called variants will be 
already normalized and filtered, ready for the downstream genotype 
analysis by the Stargazer program.

*gb* is the genome build ('hg19' or 'hg38'). *tg* is the target gene 
(e.g. 'cyp2d6') or region (e.g. 'chr22:42512500-42551883'). 
*fa* is the reference FASTA file. *bam* is the input BAM file(s). 

If you have many input BAM files, you may want to use the ``--bd`` or 
``--bl`` argument instead of manually listing individual files for *bam*.

.. note::

    BCFtools must be pre-installed.

Options
-------

--bd DIR    directory containing BAM files
--bl FILE   list of BAM files, one file per line
-o FILE     output to FILE [stdout]

genotype command
================

Synopsis
--------

pypgx genotype *[options] fa dt gb tg out [bam [bam ...]]*

Description
-----------

Call star alleles from BAM file(s).

This command runs the Stargazer genotyping pipeline without the need for 
Sun Grid Engine (SGE). It uses the ``bam2vcf`` tool to create the input 
VCF file (which is essentially a wrapper of, and therefore requires, the 
BCFtools program) and the ``bam2gdf`` tool to create the input GDF file. 
It then runs the Stargazer program to perform genotype analysis.

In order to detect strctural variation, Stargazer needs read depth data 
(i.e. a GDF file) for copy number analysis. Providing the optional 
argument ``--cg`` will generate a GDF file. If this argument is not 
used, Stargazer will run as VCF-only mode.

*fa* is the reference FASTA file. *dt* is the sequencing data type; 
use 'wgs' for whole genome sequencing data and 'ts' for targeted sequencing 
data. *gb* is the genome build ('hg19' or 'hg38'). *tg* is the target gene 
(e.g. 'cyp2d6'). *out* is the output project directory. *bam* is the 
input BAM file(s).

If you have many input BAM files, you may want to use the ``--bd`` or 
``--bl`` argument instead of manually listing individual files for *bam*.

.. note::

    BCFtools and Stargazer must be pre-installed.

Options
-------

--bd DIR    directory containing BAM files
--bl FILE   list of BAM files, one file per line
--cg STR    control gene (e.g. 'vdr') or region (e.g. 'chr12:48232319-48301814')

gt2pt command
=============

Synopsis
--------

pypgx snp *[options] gt*

Description
-----------

Call phenotypes from star alleles.

*gt* is the genotype file from Stargazer. See the API section for details.

Options
-------

-o FILE     output to FILE [stdout]

xgep command
============

Synopsis
--------

pypgx xgep *[options] conf*

Description
-----------

Run per-project genotyping for multiple genes with SGE (1).

This command runs the per-project genotyping pipeline by submitting 
jobs to the Sun Grid Engine (SGE) cluster. This is essentially an 
extension of the ``sgep`` command to genotype multiple genes.

*conf* is the configuration file. See the API section for details.

.. note::

    SGE and Stargazer must be pre-installed.

Options
-------

There are no options.