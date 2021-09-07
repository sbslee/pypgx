API
***

This section describes application programming interface (API) for PyPGx.

Archive file, semantic type, and metadata
=========================================

Each archive file (``.zip``) has a metadata file (``metadata.txt``) and a data file (e.g. ``data.tsv``, ``data.vcf``). A metadata file contains information about the data being contained within the archive, which is expressed as pairs of ``=``-separated keys and values (e.g. ``Assembly=GRCh37``):

.. list-table::
    :widths: 20 40 40
    :header-rows: 1

    * - Metadata
      - Description
      - Examples
    * - ``Assembly``
      - Reference genome assembly.
      - ``GRCh37``, ``GRCh38``
    * - ``Control``
      - Control gene.
      - ``VDR``, ``chr1:10000-20000``
    * - ``Gene``
      - Target gene.
      - ``CYP2D6``, ``GSTT1``
    * - ``Platform``
      - NGS platform.
      - ``WGS``, ``Targeted``
    * - ``Program``
      - Name of the phasing program.
      - ``Beagle``
    * - ``SemanticType``
      - Semantic type of the archive.
      - ``CovFrame[CopyNumber]``, ``Model[CNV]``

Notably, all archive files have defined semantic types, which allows us to ensure that the data that is passed to a method/command is meaningful for the operation that will be performed. Below is a list of currently defined semantic types:

- ``CovFrame[CopyNumber]``
    * CovFrame for storing target gene's per-base copy number which is computed from read depth with control statistics.
    * Requires following metadata: ``Gene``, ``Assembly``, ``SemanticType``, ``Control``.
- ``CovFrame[ReadDepth]``
    * CovFrame for storing target gene's per-base read depth which is computed from BAM files.
    * Requires following metadata: ``Gene``, ``Assembly``, ``SemanticType``.
- ``Model[CNV]``
    * Model for calling CNV in target gene.
    * Requires following metadata: ``Gene``, ``Assembly``, ``SemanticType``, ``Control``.
- ``TSV[Alleles]``
    * TSV table for storing target gene's candidate star alleles for each sample.
    * Requires following metadata: ``Gene``, ``Assembly``, ``SemanticType``, ``Program``.
- ``TSV[CNVCalls]``
    * TSV table for storing target gene's CNV call for each sample.
    * Requires following metadata: ``Gene``, ``Assembly``, ``SemanticType``, ``Control``.
- ``TSV[Statistcs]``
    * TSV table for storing control gene's various statistics on read depth. Used for converting target gene's read depth to copy number.
    * Requires following metadata: ``Control``, ``Assembly``, ``SemanticType``.
- ``VcfFrame[Consolidated]``
    * VcfFrame for storing target gene's consolidated variant data.
    * Requires following metadata: ``Gene``, ``Assembly``, ``SemanticType``, ``Program``.
- ``VcfFrame[Imported]``
    * VcfFrame for storing target gene's raw variant data.
    * Requires following metadata: ``Gene``, ``Assembly``, ``SemanticType``.
- ``VcfFrame[Phased]``
    * VcfFrame for storing target gene's phased variant data.
    * Requires following metadata: ``Gene``, ``Assembly``, ``SemanticType``, ``Program``.

pypgx.api.utils
===============

.. automodule:: pypgx.api.utils
   :members:
