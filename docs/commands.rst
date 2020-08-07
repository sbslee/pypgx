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

*snp_caller* is the SNP caller ('gatk' or 'bcftools'). *fasta_file* is the 
reference FASTA file. *target_gene* is the name of target gene (e.g. 
'cyp2d6'). *genome_build* is the genome build ('hg19' or 'hg38'). 
*data_type* is the type of sequencing data ('wgs' or 'ts'). Output files 
will be written to *proj_dir*. *bam_file* are the input BAM files.

.. warning::
    Stargazer and GATK/BCFtools must be pre-installed.

Options
-------

-h, --help          show command-specific help message and exit
--bam_dir DIR       use all BAM files in this directory as input
--bam_list FILE     list of input BAM files, one file per line
--control_gene STR  name or region of control gene (e.g. ‘vdr’,
                    ‘chr12:48232319-48301814’)
--dbsnp_file FILE   dbSNP VCF file, used by GATK to add rs numbers
--temp_dir DIR      temporary files will be written to this directory
--plot              output copy number plots.

bam2gt2 command
===============

Synopsis
--------

pypgx bam2gt2 *[options] conf_file*

Description
-----------

Convert BAM files to a genotype file [SGE].

This command runs the entire genotyping pipeline for BAM files 
with the Sun Grid Engine (SGE) cluster. By default, it will genotype 
all genes currently targeted by the Stargazer program (you can specify 
select genes too). For each gene, the command runs under the hood 
``bam2vcf`` with ``bcftools`` caller (i.e. BCFtools) or ``bam2vcf2`` 
(i.e. GATK) to create the input VCF file. The input GDF file is 
created with ``bam2gdf``.

*conf_file* is the configuration file.

.. warning::

    SGE, Stargazer and BCFtools/GATK must be pre-installed.

Options
-------

-h, --help  show command-specific help message and exit

This is what a typical configuration file for ``bam2gt2`` looks like:

    .. code-block:: python

        # File: example_conf.txt
        # To execute:
        #   $ pypgx bam2gt2 example_conf.txt
        #   $ sh ./myproject/example-qsub.sh

        # Do not make any changes to this section.
        [DEFAULT]
        control_gene = NONE
        dbsnp_file = NONE
        java_options = NONE
        plot = FALSE
        qsub_options = NONE
        sample_list = NONE
        target_genes = ALL

        # Make any necessary changes to this section.
        [USER]
        bam_list = bam-list.txt
        control_gene = vdr
        data_type = wgs
        fasta_file = hs37d5.fa
        genome_build = hg19
        project_path = ./myproject
        qsub_options = -l mem_requested=2G
        snp_caller = gatk
        target_genes = cyp2b6, cyp2d6

This table summarizes the configuration parameters specific to ``bam2gt2``:

    .. list-table::
        :widths: 25 75
        :header-rows: 1

        * - Parameter
          - Summary
        * - bam_list
          - List of input BAM files, one file per line.
        * - control_gene
          - Control gene or region.
        * - data_type
          - Data type ('wgs' or 'ts').
        * - dbsnp_file
          - dbSNP VCF file.
        * - fasta_file
          - Reference FASTA file.
        * - genome_build
          - Genome build ('hg19' or 'hg38').
        * - java_options
          - Java-specific arguments for GATK (e.g. ‘-Xmx4G’).
        * - plot
          - Output copy number plots.
        * - project_path
          - Output project directory.
        * - qsub_options
          - Options for qsub command (e.g. '-l mem_requested=2G').
        * - sample_list
          - List of samples used for inter-sample normalization 
            (e.g. 'gstt1, sample1, sample2 | ugt2b17, sample3'). 
        * - snp_caller
          - SNP caller (‘gatk’ or ‘bcftools’).
        * - target_genes
          - Names of target genes (e.g. 'cyp2d6').

gt2pt command
=============

Synopsis
--------

pypgx snp *[options] gt*

Description
-----------

Call phenotypes from star alleles.

This command is just a wrapper for the ``phenotyper`` module. See the API 
section for details.

*gt* is the genotype file. 

Options
-------

-h, --help  show command-specific help message and exit
-o FILE     output to FILE [stdout]

bam2vcf command
===============

Synopsis
--------

| pypgx bam2vcf *[options]* \\
|   *snp_caller* \\
|   *fasta_file* \\
|   *target_gene* \\
|   *output_file* \\
|   *genome_build* \\
|   *[bam_file [bam_file ...]]*

Description
-----------

Convert BAM files to a VCF file.

This command creates a single- or multi-sample VCF file from one or 
more input BAM files. The output VCF file will only contain variants 
within the target gene or region. The command is essentially a wrapper 
for the Genome Analysis Toolkit (GATK) and the BCFtools program with 
pre-specified parameters. This means the called variants will be 
already normalized and filtered, ready for the downstream genotype 
analysis by the Stargazer program.

*snp_caller* is the SNP caller ('gatk' or 'bcftools'). *fasta_file* is the 
reference FASTA file. *target_gene* is the name or region of target gene 
(e.g. 'cyp2d6', 'chr22:42512500-42551883'). Output will be written to 
*output_file*. *genome_build* is the genome build ('hg19' or 'hg38'). 
*bam_file* are the input BAM files.

.. warning::
    GATK and/or BCFtools must be pre-installed.

.. note::
    Generally, GATK is more accurate but much slower than BCFtools. 
    For instance, SNP calling for 70 WGS samples for the CYP2D6 gene 
    takes 19 min with the ``gatk`` caller but only 2 min with the 
    ``bcftools`` caller. Therefore, if you have many samples and you do 
    not have access to Sun Grid Engine (SGE) for parallelism, we 
    recommend that you use ``bcftools``. If you have SGE and want to 
    use GATK, please check ``bam2vcf2``.

Options
-------

-h, --help         show command-specific help message and exit
--bam_dir DIR      use all BAM files in this directory as input
--bam_list FILE    list of input BAM files, one file per line
--dbsnp_file FILE  dbSNP VCF file, used by GATK to add rs numbers
--java_options STR  Java-specific arguments for GATK (e.g. '-Xmx4G')
--temp_dir DIR     temporary files will be written to this directory

bam2vcf2 command
================

Synopsis
--------

pypgx bam2vcf2 *[options] conf_file*

Description
-----------

Convert BAM files to a VCF file [SGE].

This command outputs a single- or multi-sample VCF file from one or 
more input BAM files. The output VCF file will only contain variants 
within the target gene or region. This command is essentially a 
wrapper with pre-specified parameters for the Genome Analysis Toolkit 
(GATK). It also uses Sun Grid Engine (SGE) for parallelism to make 
GATK run faster.

*conf_file* is the configuration file.

.. warning::
    GATK and SGE must be pre-installed.

Options
-------

-h, --help  show command-specific help message and exit

This is what a typical configuration file for ``bam2vcf2`` looks like:

    .. code-block:: python

        # File: example_conf.txt
        # To execute:
        #   $ pypgx bam2vcf2 example_conf.txt
        #   $ sh ./myproject/example-qsub.sh

        # Do not make any changes to this section.
        [DEFAULT]
        dbsnp_file = NONE
        java_options = NONE
        qsub_options = NONE

        # Make any necessary changes to this section.
        [USER]
        bam_list = bam-list.txt
        dbsnp_file = dbsnp.vcf
        fasta_file = reference.fa
        genome_build = hg19
        java_options = -Xmx4G
        project_path = ./myproject
        qsub_options = -l mem_requested=4G
        target_gene = cyp2d6

This table summarizes the configuration parameters specific to ``bam2vcf2``:

    .. list-table::
       :widths: 25 75
       :header-rows: 1

       * - Parameter
         - Summary
       * - bam_list
         - List of input BAM files, one file per line.
       * - dbsnp_file
         - dbSNP VCF file.
       * - fasta_file
         - Reference FASTA file.
       * - genome_build
         - Genome build ('hg19' or 'hg38').
       * - java_options
         - Java-specific arguments for GATK (e.g. ‘-Xmx4G’).
       * - project_path
         - Output project directory.
       * - qsub_options
         - Options for qsub command (e.g. '-l mem_requested=2G').
       * - target_gene
         - Name of target gene (e.g. 'cyp2d6'). 
           Also accepts a BED file.

bam2gdf command
===============

Synopsis
--------

| pypgx bam2gdf *[options]* \\
|   *genome_build* \\
|   *target_gene* \\
|   *control_gene* \\
|   *output_file* \\
|   *[bam_file [bam_file ...]]*

Description
-----------

Convert BAM files to a GDF file.

This command calculates read depth from BAM files and then outputs a
GDF (GATK-DepthOfCoverage Format) file, which is one of the input 
files for the Stargazer program. Even though ``gatk DepthOfCoverage`` 
could still be used to make GDF files, we recommend that you use this 
command because the former is too heavy (i.e. requires too much memory) 
for such a simple task (i.e. counting reads). The latter uses 
``samtools depth`` under the hood, which is way faster and requires 
way less memory. Another nice about using ``bam2gdf`` instead of 
``samtools depth`` is that everything is already parametrized for 
compatibility with Stargazer. 

*genome_build* is the genome build ('hg19' or 'hg38'). *target_gene* is 
the name of target gene (e.g. 'cyp2d6'). *control_gene* is the name or 
region of control gene (e.g. 'vdr', 'chr12:48232319-48301814'). Output will 
be written to *output_file*. *bam_file* are the input BAM files.

.. note::
    You do NOT need to install ``samtools`` to run this command.

Options
-------

-h, --help       show command-specific help message and exit
--bam_dir DIR    use all BAM files in this directory as input
--bam_list FILE  list of input BAM files, one file per line

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

-h, --help  show command-specific help message and exit
-o FILE     output to FILE [stdout]

bam2html command
================

Synopsis
--------

pypgx bam2html *[options] conf_file*

Description
-----------

Run per-sample genotyping for multiple genes with SGE.

This command runs the per-sample genotyping pipeline by submitting 
jobs to the Sun Grid Engine (SGE) cluster. This essentially deploys 
the ``genotype`` command to multiple genes in parallel. After genotype 
analysis is complete, it will merge the genotype results and then 
generate a HTML report using the ``gt2html`` command.

*conf_file* is the configuration file.

.. note::

    BCFtools, SGE and Stargazer must be pre-installed.

Options
-------

-h, --help  show command-specific help message and exit

This is what a typical configuration file for ``sges`` looks like:

    .. code-block:: python

        # File: example_conf.txt
        # To execute:
        #   $ pypgx sges example_conf.txt
        #   $ sh ./myproject/example-qsub.sh

        # Do not make any changes to this section.
        [DEFAULT]
        target_genes = ALL
        control_gene = NONE
        plot = FALSE
        qsub_options = NONE

        # Make any necessary changes to this section.
        [USER]
        snp_caller = gatk
        fasta_file = reference.fa
        project_path = ./myproject
        genome_build = hg19
        data_type = wgs
        bam_file = in.bam
        qsub_options = -l mem_requested=2G
        target_genes = cyp2b6, cyp2d6
        control_gene = vdr

This table summarizes the configuration parameters specific to ``sges``:

    .. list-table::
       :widths: 25 75
       :header-rows: 1

       * - Parameter
         - Summary
       * - bam_file
         - BAM file.
       * - control_gene
         - Name or region of control gene 
           (e.g. 'vdr', 'chr12:48232319-48301814').
       * - data_type
         - Data type ('wgs' or 'ts').
       * - fasta_file
         - Reference FASTA file.
       * - genome_build
         - Genome build ('hg19' or 'hg38').
       * - plot
         - Output copy number plots.
       * - project_path
         - Output project directory.
       * - qsub_options
         - Options for qsub command (e.g. '-l mem_requested=2G').
       * - target_genes
         - Names of target genes (e.g. 'cyp2d6').

fq2bam command
==============

Synopsis
--------

pypgx fq2bam *[options] conf_file*

Description
-----------

Create BAM file(s) from FASTQ file(s).

*conf_file* is the configuration file.

Options
-------

-h, --help  show command-specific help message and exit

This is what a typical configuration file for ``fq2bam`` looks like:

    .. code-block:: python

        # File: example_conf.txt
        # Do not make any changes to this section.
        [DEFAULT]
        platform = illumina
        qsub_options1 = NONE
        qsub_options2 = NONE
        read_length = 150
        threads = 1

        # Make any necessary changes to this section.
        [USER]
        bed_file = in.bed
        fasta_file = reference.fa
        library = awesome_experiment
        manifest_file = manifest.txt
        project_path = /path/to/project/
        qsub_options1 = -V -q biall.q -S /bin/bash -pe pePAC 15
        qsub_options2 = -V -q biall.q -S /bin/bash
        threads = 15
        vcf_files = in1.vcf, in2.vcf, in3.vcf

This table summarizes the configuration parameters specific to ``fq2bam``:

    .. list-table::
        :widths: 25 75
        :header-rows: 1

        * - Parameter
          - Summary
        * - bed_file
          - BED file.
        * - fasta_file
          - Reference FASTA file.
        * - library
          - Sequencing library name.
        * - manifest_file
          - Manifest file.
        * - platform
          - Sequencing platform.
        * - project_path
          - Output project directory.
        * - qsub_options1
          - Options for the first qsub command. Recommended to set a parallel environment.
        * - qsub_options2
          - Options for the second qsub command.
        * - read_length
          - Sequence read length.
        * - threads
          - Number of threads.
        * - vcf_files
          - Reference VCF files used for base quality score recalibration.

bam2bam command
===============

Synopsis
--------

pypgx bam2bam *[options] conf_file*

Description
-----------

Realign BAM files to another reference genome [SGE].

*conf_file* is the configuration file.

Options
-------

-h, --help  show command-specific help message and exit

This is what a typical configuration file for ``bam2bam`` looks like:

    .. code-block:: python

        # File: example_conf.txt
        # Do not make any changes to this section.
        [DEFAULT]
        java_heap = -Xmx2g
        platform = illumina
        qsub_options1 = NONE
        qsub_options2 = NONE
        threads = 1

        # Make any necessary changes to this section.
        [USER]
        fasta_file = reference.fa
        gatk_tool = GenomeAnalysisTK.jar
        library = awesome_experiment
        manifest_file = manifest.txt
        picard_tool = picard.jar
        project_path = /path/to/project/
        qsub_options1 = -q nick-grad.q -l mem_requested=2G -pe serial 1
        qsub_options2 = -q nick-grad.q -l mem_requested=2G
        vcf_files = in1.vcf, in2.vcf, in3.vcf

This table summarizes the configuration parameters specific to ``bam2bam``:

    .. list-table::
        :widths: 25 75
        :header-rows: 1

        * - Parameter
          - Summary
        * - fasta_file
          - Reference FASTA file.
        * - gatk_tool
          - GATK program.
        * - java_heap
          - Java heap size.
        * - library
          - Sequencing library name.
        * - manifest_file
          - Manifest file.
        * - picard_tool
          - Picard program.
        * - platform
          - Sequencing platform.
        * - project_path
          - Output project directory.
        * - qsub_options1
          - Options for the first qsub command. Recommended to set a parallel environment.
        * - qsub_options2
          - Options for the second qsub command.
        * - threads
          - Number of threads.
        * - vcf_files
          - Reference VCF files used for base quality score recalibration.

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

-h, --help  show command-specific help message and exit
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

-h, --help  show command-specific help message and exit
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

-h, --help   show command-specific help message and exit
-o FILE      output to FILE [stdout]
--test_mode  only extract first three guidelines for testing

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

-h, --help  show command-specific help message and exit
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

-h, --help  show command-specific help message and exit
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

-h, --help  show command-specific help message and exit
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

-h, --help  show command-specific help message and exit
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

-h, --help  show command-specific help message and exit
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

-h, --help  show command-specific help message and exit
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

-h, --help  show command-specific help message and exit

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

-h, --help  show command-specific help message and exit

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

-h, --help  show command-specific help message and exit
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

-h, --help  show command-specific help message and exit
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

-h, --help  show command-specific help message and exit
-o FILE     output to FILE [stdout]

compare2 command
================

Synopsis
--------

pypgx compare2 *[options] truth_file test_file sample_map*

Description
-----------

Compare two genotype files.

*truth_file* is the truth genotype file. *test_file* is the test genotype 
file. *sample_map* is the tab-delimited text file for sample name mapping.

Options
-------

-h, --help  show command-specific help message and exit
-o FILE     output to FILE [stdout]

compvcf command
===============

Synopsis
--------

pypgx compvcf *[options] truth_file test_file sample_map*

Description
-----------

Compare two VCF files.

*truth_file* is the truth VCF file. *test_file* is the test VCF 
file. *sample_map* is the tab-delimited text file for sample name mapping.

Options
-------

-h, --help  show command-specific help message and exit
