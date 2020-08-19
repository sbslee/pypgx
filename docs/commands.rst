Commands
********

This section describes how to use PyPGx as a command-line program.
For the public API of Python module ``pypgx``, please see the API section.

bam2gt command
==============

Convert BAM files to a genotype file.

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

Positional arguments
--------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *snp_caller*
     - SNP caller ('gatk' or 'bcftools').
   * - *fasta_file*
     - Reference FASTA file.
   * - *target_gene*
     - Name of target gene (e.g. 'cyp2d6').
   * - *genome_build*
     - Genome build ('hg19' or 'hg38').
   * - *data_type*
     - Type of sequencing data ('wgs' or 'ts').
   * - *proj_dir*
     - Output files will be written to *proj_dir*.
   * - *bam_file*
     - Input BAM files.

Optional arguments
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *-h, \--help*
     - Show command-specific help message and exit.
   * - *\\-\\-bam_dir DIR*
     - Use all BAM files in this directory as input.
   * - *\--bam_list FILE*
     - List of input BAM files, one file per line.
   * - *\--control_gene STR*
     - Name or region of control gene (e.g. ‘vdr’, ‘chr12:48232319-48301814’).
   * - *\--dbsnp_file FILE*
     - dbSNP VCF file, used by GATK to add rs numbers.
   * - *\--temp_dir DIR*
     - Temporary files will be written to this directory.
   * - *\--plot*
     - Output copy number plots.

Description
-----------

This command runs the entire genotyping pipeline for BAM files, 
without the need for Sun Grid Engine (SGE). Under the hood, it 
uses the ``bam2vcf`` command to create the input VCF file and 
the ``bam2gdf`` command to create the input GDF file. It then 
performs genotype analysis using the Stargazer program.

In order to detect strctural variation, Stargazer requires read 
depth data (i.e. a GDF file) for copy number analysis. Providing 
the optional argument ``--control_gene`` will generate a GDF file. 
If this argument is not provided, Stargazer will run as VCF-only mode.

.. warning::
    Stargazer and GATK/BCFtools must be pre-installed.

bam2gt2 command
===============

Convert BAM files to a genotype file [SGE].

Synopsis
--------

pypgx bam2gt2 *[options] conf_file*

Positional arguments
--------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *conf_file*
     - Configuration file.

Optional arguments
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *-h, \--help*
     - Show command-specific help message and exit.

Description
-----------

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

Call phenotypes from star alleles.

Synopsis
--------

pypgx gt2pt *[options] gt*

Positional arguments
--------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *gt*
     - Genotype file.

Optional arguments
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *-h, \--help*
     - Show command-specific help message and exit.
   * - *\-o*
     - Output to FILE [stdout].

Description
-----------

This command is just a wrapper for the ``phenotyper`` module. See the API 
section for details.

bam2vcf command
===============

Convert BAM files to a VCF file.

Synopsis
--------

| pypgx bam2vcf *[options]* \\
|   *snp_caller* \\
|   *fasta_file* \\
|   *target_gene* \\
|   *output_file* \\
|   *genome_build* \\
|   *[bam_file [bam_file ...]]*

Positional arguments
--------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *snp_caller*
     - SNP caller ('gatk' or 'bcftools').
   * - *fasta_file*
     - Reference FASTA file.
   * - *target_gene*
     - Name or region of target gene (e.g. 'cyp2d6', 'chr22:42512500-42551883').
   * - *output_file*
     - Output will be written to *output_file*.
   * - *genome_build*
     - Genome build ('hg19' or 'hg38').
   * - *bam_file*
     - Input BAM files.

Optional arguments
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *-h, \--help*
     - Show command-specific help message and exit.
   * - *\--bam_dir DIR*
     - Use all BAM files in this directory as input.
   * - *\--bam_list FILE*
     - List of input BAM files, one file per line.
   * - *\--dbsnp_file FILE*
     - dbSNP VCF file, used by GATK to add rs numbers.
   * - *\--java_options STR*
     - Java-specific arguments for GATK (e.g. '-Xmx4G').
   * - *\--temp_dir DIR*
     - Temporary files will be written to this directory.

Description
-----------

This command creates a single- or multi-sample VCF file from one or 
more input BAM files. The output VCF file will only contain variants 
within the target gene or region. The command is essentially a wrapper 
for the Genome Analysis Toolkit (GATK) and the BCFtools program with 
pre-specified parameters. This means the called variants will be 
already normalized and filtered, ready for the downstream genotype 
analysis by the Stargazer program.


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

bam2vcf2 command
================

Convert BAM files to a VCF file [SGE]

Synopsis
--------

pypgx bam2vcf2 *[options] conf_file*

Positional arguments
--------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *conf_file*
     - Configuration file.

Optional arguments
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *-h, \--help*
     - Show command-specific help message and exit.

Description
-----------

This command outputs a single- or multi-sample VCF file from one or 
more input BAM files. The output VCF file will only contain variants 
within the target gene or region. This command is essentially a 
wrapper with pre-specified parameters for the Genome Analysis Toolkit 
(GATK). It also uses Sun Grid Engine (SGE) for parallelism to make 
GATK run faster.

*conf_file* is the configuration file.

.. warning::
    GATK and SGE must be pre-installed.

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

Convert BAM files to a GDF file.

Synopsis
--------

| pypgx bam2gdf *[options]* \\
|   *genome_build* \\
|   *target_gene* \\
|   *control_gene* \\
|   *output_file* \\
|   *[bam_file [bam_file ...]]*

Positional arguments
--------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *genome_build*
     - Genome build ('hg19' or 'hg38').
   * - *target_gene*
     - Name of target gene (e.g. 'cyp2d6').
   * - *control_gene*
     - Name or region of control gene (e.g. 'vdr', 'chr12:48232319-48301814').
   * - *output_file*
     - Output will be written to *output_file*.
   * - *bam_file*
     - Input BAM files.

Optional arguments
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *-h, \--help*
     - Show command-specific help message and exit.
   * - *\--bam_dir DIR*
     - Use all BAM files in this directory as input.
   * - *\--bam_list FILE*
     - List of input BAM files, one file per line.

Description
-----------

This command converts BAM files to a GDF file.

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

.. note::
    You do NOT need to install ``samtools`` to run this command.

gt2html command
===============

Create HTML report using Stargazer data.

Synopsis
--------

pypgx gt2html *[options] gt*


Positional arguments
--------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *gb*
     - Genotype file.

Optional arguments
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *-h, \--help*
     - Show command-specific help message and exit.
   * - *-o*
     - Output to FILE [stdout].

Description
-----------

This command creates HTML report using Stargazer data.

bam2html command
================

Run per-sample genotyping for multiple genes with SGE.

Synopsis
--------

pypgx bam2html *[options] conf_file*

Positional arguments
--------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *conf_file*
     - Configuration file.

Optional arguments
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *-h, \--help*
     - Show command-specific help message and exit.

Description
-----------

This command runs per-sample genotyping for multiple genes with SGE.

This command runs the per-sample genotyping pipeline by submitting 
jobs to the Sun Grid Engine (SGE) cluster. This essentially deploys 
the ``genotype`` command to multiple genes in parallel. After genotype 
analysis is complete, it will merge the genotype results and then 
generate a HTML report using the ``gt2html`` command.

.. note::

    BCFtools, SGE and Stargazer must be pre-installed.

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

Create BAM file(s) from FASTQ file(s).

Synopsis
--------

pypgx fq2bam *[options] conf_file*

Positional arguments
--------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *conf_file*
     - Configuration file.

Optional arguments
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *-h, \--help*
     - Show command-specific help message and exit.

Description
-----------

This command creates BAM file(s) from FASTQ file(s).

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

Realign BAM files to another reference genome [SGE].

Synopsis
--------

pypgx bam2bam *[options] conf_file*

Positional arguments
--------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *conf_file*
     - Configuration file.

Optional arguments
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *-h, \--help*
     - Show command-specific help message and exit.

Description
-----------

This command realign BAM files to another reference genome using SGE.

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

Create SDF file from BAM file(s).

Synopsis
--------

pypgx bam2sdf *[options] gb tg cg bam [bam ...]*

Positional arguments
--------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *gb*
     - Genome build ('hg19' or 'hg38').
   * - *tg*
     - Target gene (e.g. 'cyp2d6').
   * - *cg*
     - Control gene (e.g. 'vdr') or region (e.g. 'chr12:48232319-48301814').
   * - *bam*
     - BAM file.

Optional arguments
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *-h, \--help*
     - Show command-specific help message and exit.
   * - *-o*
     - Output to FILE [stdout].

Description
-----------

This command creates SDF file from BAM files.

sdf2gdf command
===============

Create GDF file from SDF file.

Synopsis
--------

pypgx sdf2gdf *[options] sdf id [id ...]*

Positional arguments
--------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *sdf*
     - SDF file.
   * - *id*
     - Sample ID.

Optional arguments
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *-h, \--help*
     - Show command-specific help message and exit.
   * - *-o*
     - Output to FILE [stdout].

Description
-----------

This command creates GDF file from SDF file.

pgkb command
============

Extract CPIC guidelines using PharmGKB API.

Synopsis
--------

pypgx pgkb *[options]*

Positional arguments
--------------------

There are no positional arguments for this command.

Optional arguments
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *\--test_mode*
     - Only extract first three guidelines for testing.

Description
-----------

This command extracts CPIC recommendations for prescription drugs using 
PharmGKB API.

minivcf command
===============

Slice VCF file.

Synopsis
--------

pypgx minivcf *[options] vcf_file region*

Positional arguments
--------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *vcf_file*
     - VCF file.
   * - *region*
     - Target region.

Optional arguments
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *-h, \--help*
     - Show command-specific help message and exit.

Description
-----------

This command slices a VCF file for the given region.

mergevcf command
================

Merge VCF files.

Synopsis
--------

pypgx mergevcf *[options] vcf_file [vcf_file ...]*

Positional arguments
--------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *vcf_file*
     - VCF files to be merged.

Optional arguments
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *-h, \--help*
     - Show command-specific help message and exit.
   * - *\--region*
     - Target region.

Description
-----------

This command merges VCF files with single sample. It's assumed that the VCF 
files share the same variant sites. In the upcoming version, these 
restrictions will be lifted and the command will be able to merge VCF files 
with any number of samples and with different sets of variants.

summary command
===============

Create summary file using Stargazer data.

Synopsis
--------

pypgx summary *[options] gt*

Positional arguments
--------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *gt*
     - Genotype file.

Optional arguments
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *-h, \--help*
     - Show command-specific help message and exit.
   * - *-o*
     - Output to FILE [stdout].

Description
-----------

This command creates summary file using Stargazer data.

meta command
============

Create meta file from summary files.

Synopsis
--------

pypgx meta *[options] sf [sf ...]*

Positional arguments
--------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *sf*
     - Summary file from the ``summary`` command.


Optional arguments
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *-h, \--help*
     - Show command-specific help message and exit.
   * - *-o*
     - Output to FILE [stdout].

Description
-----------

This command creates meta comparison file from summary files.

compare command
===============

Compare genotype files.

Synopsis
--------

pypgx compare *[options] gt [gt ...]*

Positional arguments
--------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *gt*
     - Genotype file.


Optional arguments
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *-h, \--help*
     - Show command-specific help message and exit.
   * - *-o*
     - Output to FILE [stdout].

Description
-----------

This command can compare multiple genotype files at once.

check command
=============

Checks table files for Stargazer.

Synopsis
--------

| pypgx check *[options]* \\
|   *star_table* \\
|   *snp_table* \\

Positional arguments
--------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *star_table*
     - Star allele table file (``star_table.txt``).
   * - *snp_table*
     - SNP table file (``snp_table.txt``).

Optional arguments
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *-h, \--help*
     - Show command-specific help message and exit.

Description
-----------

This command is meant to be used for Stargazer development.

liftover command
================

Convert variants in SNP table from hg19 to hg38.

Synopsis
--------

| pypgx liftover *[options]* \\
|   *star_table* \\
|   *snp_table* \\
|   *target_gene*

Positional arguments
--------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *star_table*
     - Star allele table file (``star_table.txt``).
   * - *snp_table*
     - SNP table file (``snp_table.txt``).
   * - *target_gene*
     - Target gene.

Optional arguments
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *-h, \--help*
     - Show command-specific help message and exit.

Description
-----------

This command is meant to be used for Stargazer development.

peek command
============

Find all possible star alleles from VCF file.

Synopsis
--------

pypgx peek *[options] vcf_file*

Positional arguments
--------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *vcf_file*
     - Stargazer VCF file (``finalized.vcf``).

Optional arguments
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *-h, \--help*
     - Show command-specific help message and exit.

Description
-----------

This command returns summary of the status of all possibile star alleles 
that can be called from the VCF file.

viewsnp command
===============

View SNP data for pairs of sample/star allele.

Synopsis
--------

| pypgx viewsnp *[options]* \\
|   *vcf_file* \\
|   *query [query ...]*

Positional arguments
--------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *vcf_file*
     - Stargazer VCF file (``finalized.vcf``).
   * - *query [query ...]*
     - Pair of sample and star allele separated by ``/`` 
       (e.g. ``SAMPLE1/*4``).

Optional arguments
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *-h, -\-help*
     - Show command-specific help message and exit.

Description
-----------

This command shows the SNP data for given pairs of a sample and a star 
allele. It's designed to be used after running Stargazer.

Here's a complete example with real NGS data.

.. code-block:: python

   # Install Stargazer.
   python -m pip install git+https://github.com/sbslee/stargazer

   # Download example data.
   git clone https://github.com/sbslee/stargazer
   cd stargazer/example

   # Run Stargazer as in:
   # https://stargazer.readthedocs.io/en/latest/tutorial.html#example-1.
   stargazer \
     wgs \
     hg19 \
     cyp2d6 \
     getrm-cyp2d6-vdr.joint.filtered.vcf \
     ./ex1-getrm-cyp2d6-vdr \
     --gdf getrm-cyp2d6-vdr.gdf \
     --cg vdr

   # Run viewsnp.
   pypgx viewsnp \
     ex1-getrm-cyp2d6-vdr/finalized.vcf \
     316ab006177d41b484982d7fa4d851ad/*21 \
     2c9f234af49b4f6a970d8ddef07358e5/*4

The output will look like this::

    <sample=316ab006177d41b484982d7fa4d851ad,star=*21>
    hg19_pos	wt_allele	var_allele	hg19_allele	type	so	impact	effect	hap1_allele	hap2_allele	gt	hap1_ad	hap2_ad	hap1_af	hap2_af
    42522613	C	G	G	tag	missense_variant	low_impact	S486T	C	G	0|1	19	10	0.66	0.34
    42523409	T	G	G	tag	intron_variant	low_impact	no_effect	T	G	0|1	19	23	0.45	0.55
    42523943	G	A	A	tag	missense_variant	low_impact	R296C	G	A	0|1	21	15	0.58	0.42
    42524213	C	CG	C	core	frameshift_variant	high_impact	frameshift	C	CG	0|1	14	12	0.54	0.46
    42525132	C	G	G	tag	synonymous_variant	low_impact	V136#	C	G	0|1	18	28	0.39	0.61
    42526580	C	G	G	tag	intron_variant	low_impact	no_effect	C	G	0|1	22	23	0.49	0.51
    42528382	G	C	C	tag	upstream_gene_variant	low_impact	no_effect	G	C	0|1	14	14	0.50	0.50
    <sample=2c9f234af49b4f6a970d8ddef07358e5,star=*4>
    hg19_pos	wt_allele	var_allele	hg19_allele	type	so	impact	effect	hap1_allele	hap2_allele	gt	hap1_ad	hap2_ad	hap1_af	hap2_af
    42524947	C	T	C	core	splice_acceptor_variant	high_impact	splicing_defect	T	C	1|0	14	23	0.38	0.62
    42526694	G	A	G	tag	missense_variant	high_impact	P34S	A	G	1|0	26	16	0.62	0.38

compgt command
==============

Compute the concordance between two genotype files.

Synopsis
--------

| pypgx compgt *[options]* \\
|   *truth_file* \\
|   *test_file* \\
|   *sample_map*

Positional arguments
--------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *truth_file*
     - Truth genotype file.
   * - *test_file*
     - Test genotype file.
   * - *sample_map*
     - Tab-delimited text file with two columns representing 
       the truth and test sample names.

Optional arguments
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *-h, \--help*
     - Show command-specific help message and exit.

Description
-----------

This command computes the concordance between genotype data (e.g. 
``*1/*4``) of one samples in each of the genotype files, one being 
considered the truth and the other being the test.

compvcf command
===============

Calculate the concordance between two VCF files.

Synopsis
--------

| pypgx compvcf *[options]* \\
|   *truth_file* \\
|   *test_file* \\
|   *sample_map*

Positional arguments
--------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *truth_file*
     - Truth VCF file.
   * - *test_file*
     - Test VCF file.
   * - *sample_map*
     - Tab-delimited text file with two columns representing 
       the truth and test sample names.

Optional arguments
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *-h, \--help*
     - Show command-specific help message and exit.

Description
-----------

This command calculates the concordance between genotype data (e.g. ``0/1``) 
of one samples in each of the VCF files, one being considered the truth and 
the other being the test. The concordance is broken into separate results 
sections for SNP and Indel. Summary and detailed statistics are reported.

Please note that the comparison is restricted to sites that are biallelic and 
have no missing genotypes (e.g. ``./.``).

This table summarizes the column headers of the output.

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Header
     - Summary
   * - name1
     - Truth sample name.
   * - name2
     - Test sample name.
   * - snv_tn
     - Number of true negatives for SNV.
   * - snv_tp
     - Number of true positives for SNV.
   * - snv_fn
     - Number of false negatives for SNV.
   * - snv_fp
     - Number of false positives for SNV.
   * - snv_tpr
     - True positive rate for SNV.
   * - snv_tnr
     - True negative rate for SNV.
   * - snv_con
     - Accuracy for SNV.
   * - indel_tn
     - Number of true negatives for Indel.
   * - indel_tp
     - Number of true positives for Indel.
   * - indel_fn
     - Number of false negatives for Indel.
   * - indel_fp
     - Number of false positives for Indel.
   * - indel_tpr
     - True positive rate for Indel.
   * - indel_tnr
     - True negative rate for Indel.
   * - indel_con
     - Accuracy for Indel.
   * - all_tn
     - Number of true negatives for SNV+Indel.
   * - all_tp
     - Number of true positives for SNV+Indel.
   * - all_fn
     - Number of false negatives for SNV+Indel.
   * - all_fp
     - Number of false positives for SNV+Indel.
   * - all_tpr
     - True positive rate for SNV+Indel.
   * - all_tnr
     - True negative rate for SNV+Indel.
   * - all_con
     - Accuracy for SNV+Indel.

unicov command
==============

Compute the uniformity of sequencing coverage.

Synopsis
--------

| pypgx unicov *[options]* \\
|   *bed_file* \\
|   *[bam_file [bam_file ...]]*

Positional arguments
--------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *bed_file*
     - BED file.
   * - *bam_file*
     - Input BAM files.

Optional arguments
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Summary
   * - *-h, \--help*
     - Show command-specific help message and exit.
   * - *\--bam_dir DIR*
     - Use all BAM files in this directory as input.
   * - *\--bam_list FILE*
     - List of input BAM files, one file per line.

Description
-----------

This command evaluates the uniformity of sequencing coverage by computing 
% of base pairs that were sequenced at various coverages. Only regions 
specified in the BED file are computed.
