API
***

This section describes the public API of Python module ``pypgx``. For how to 
use PyPGx as a command-line program, please see the Commands section.

Major modules
=============

pypgx.pgkb module
-----------------

.. automodule:: pypgx.pgkb
    :members:

pypgx.report module
-------------------

.. automodule:: pypgx.report
    :members:

pypgx.sdf2gdf module
--------------------

.. automodule:: pypgx.sdf2gdf
    :members:

pypgx.bam2sdf module
--------------------

.. automodule:: pypgx.bam2sdf
    :members:

pypgx.bam2gdf module
--------------------

.. automodule:: pypgx.bam2gdf
    :members:

pypgx.minivcf module
--------------------

.. automodule:: pypgx.minivcf
    :members:

pypgx.merge module
------------------

.. automodule:: pypgx.merge
    :members:

pypgx.summary module
--------------------

.. automodule:: pypgx.summary
    :members:

pypgx.meta module
-----------------

.. automodule:: pypgx.meta
    :members:

pypgx.compare module
--------------------

.. automodule:: pypgx.compare
    :members:

pypgx.remap module
------------------

.. automodule:: pypgx.remap
    :members:

pypgx.fq2bam module
-------------------

.. automodule:: pypgx.fq2bam
    :members:

pypgx.sges module
-----------------

.. automodule:: pypgx.sges
    :members:

pypgx.sgep module
-----------------

.. automodule:: pypgx.sgep
    :members:

pypgx.sgea module
-----------------

.. automodule:: pypgx.sgea
    :members:

pypgx.cpa module
----------------

.. automodule:: pypgx.cpa
    :members:

pypgx.plotcov module
--------------------

.. automodule:: pypgx.plotcov
    :members:

pypgx.check module
------------------

.. automodule:: pypgx.check
    :members:

pypgx.liftover module
---------------------

.. automodule:: pypgx.liftover
    :members:

pypgx.peek module
-----------------

.. automodule:: pypgx.peek
    :members:

pypgx.snp module
----------------

.. automodule:: pypgx.snp
    :members:

Auxiliary modules
=================

pypgx.common module
-------------------

.. automodule:: pypgx.common
    :members:
    :undoc-members:

pypgx.sglib module
------------------

.. automodule:: pypgx.sglib
    :members:
    :undoc-members:

Configuration arguments
=======================

A number of PyPGx commands (e.g. sgep) accept a configuration file as input. 
This table summarizes the configuration arguments.

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Argument
     - Summary
   * - bam_file
     - Bam file.
   * - control_gene
     - Control gene for Stargazer. [NONE]
   * - data_type
     - Input data type (wgs, ts, chip).
   * - dbsnp_file
     - dbSNP VCF file.
   * - fasta_file
     - Reference sequence file.
   * - gatk_tool
     - Path to GATK file.
   * - genome_build
     - Genome build (hg19, hg38). [hg19]
   * - manifest_file
     - Manifest file.
   * - mapping_quality
     - Minimum mapping quality used for counting reads. [1]
   * - output_prefix
     - Output prefix. [pypgx]
   * - project_path
     - Path to output project directory.
   * - qsub_options
     - Options for qsub (e.g. -V -l mem_requested=10G).
   * - stargazer_tool
     - Path to Stargazer directory.
   * - target_gene
     - Target gene for Stargazer.
   * - target_genes
     - Target genes for Stargazer. [ALL]
   * - vcf_only
     - If true, do not use read depth data. [FALSE]
