Commands
********

This section describes how to use PyPGx as a command-line program.
For the public API of Python module ``pypgx``, please see the API section.

bam2gt command
==============

Synopsis
--------

| pypgx bam2gt *[options]* \\
|   *snp_caller* \\
|   *fasta_file* \\
|   *target_gene* \\
|   *genome_build* \\
|   *data_type* \\
|   *proj_dir* \\
|   *[bam_file [bam_file ...]]*

Description
-----------

Convert BAM files to a genotype file.

This command runs the entire genotyping pipeline for BAM files, 
without the need for Sun Grid Engine (SGE). Under the hood, it 
uses the ``bam2vcf`` command to create the input VCF file and 
the ``bam2gdf`` command to create the input GDF file. It then 
performs genotype analysis using the Stargazer program.

In order to detect strctural variation, Stargazer requires read 
depth data (i.e. a GDF file) for copy number analysis. Providing 
the optional argument ``--control_gene`` will generate a GDF file. 
If this argument is not provided, Stargazer will run as VCF-only mode.

*snp_caller* is the SNP caller ('gatk' or 'bcftools'). 
*fasta_file* is the reference FASTA file.
*target_gene* is the target gene (e.g. 'cyp2d6') or 
region (e.g. ‘chr22:42512500-42551883’). 
*genome_build* is the genome build ('hg19' or 'hg38'). 
*data_type* is the input data type ('wgs' or 'ts'). 
*proj_dir* is the output project directory. 
*bam_file* is the input BAM file.

Options
-------

--bam_dir DIR       any BAM files in DIR will be used as input [null]
--bam_list FILE     list of BAM files, one file per line [null]
--control_gene STR  control gene
--dbsnp_file FILE   dbSNP VCF file
--temp_dir DIR      temporary files will be written to DIR [/tmp]

| **--bam_dir DIR**
|     any BAM files in DIR will be used as input [null]
| --bam_list FILE
|     list of BAM files, one file per line [null]
| --control_gene STR
|     control gene
| --dbsnp_file FILE
|     dbSNP VCF file
| --temp_dir DIR
|     temporary files will be written to DIR [/tmp]

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


gt2html command
===============

Synopsis
--------

pypgx gt2html *[options] gt*

Description
-----------

Create HTML report using Stargazer data.

*gt* is the genotype file from Stargazer.

Options
-------

-o FILE     output to FILE [stdout]

bam2html command
================

Synopsis
--------

pypgx bam2html *[options] conf*

Description
-----------

Run per-sample genotyping for multiple genes with SGE.

This command runs the per-sample genotyping pipeline by submitting 
jobs to the Sun Grid Engine (SGE) cluster. This essentially deploys 
the ``genotype`` command to multiple genes in parallel. After genotype 
analysis is complete, it will merge the genotype results and then 
generate a HTML report using the ``gt2html`` command.

*conf* is the configuration file. See the API section for details.

.. note::

    BCFtools, SGE and Stargazer must be pre-installed.

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

    BCFtools, SGE and Stargazer must be pre-installed.

Options
-------

There are no options.

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

    BCFtools, SGE and Stargazer must be pre-installed.

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

bam2bam command
===============

Synopsis
--------

pypgx bam2bam *[options] conf*

Description
-----------

Remap BAM file(s) to different reference.

*conf* is the configuration file. See the API section for details.

Options
-------

There are no options.

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
