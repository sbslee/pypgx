API
***

This section describes application programming interface (API) for PyPGx.

Archive and semantic type
=========================

- ``CovFrame[CopyNumber]``
    * CovFrame for storing target gene's per-base copy number which is computed from read depth with control statistics.
    * Requires following metadata: ``Gene``, ``Assembly``, ``SemanticType``, ``Control``.
- ``CovFrame[ReadDepth]``
    * CovFrame for storing target gene's per-base read depth which is computed from BAM files.
    * Requires following metadata: ``Gene``, ``Assembly``, ``SemanticType``.
- ``Model[CNV]``
    * Model for calling CNV in target gene.
    * Requires following metadata:
- ``TSV[Alleles]``
    * TSV table for storing target gene's candidate star alleles for each sample.
    * Requires following metadata:
- ``TSV[CNVCalls]``
    * TSV table for storing target gene's CNV call for each sample.
    * Requires following metadata: ``Gene``, ``Assembly``, ``SemanticType``, ``Control``.
- ``TSV[Statistcs]``
    * TSV table for storing control gene's various statistics on read depth.
    * Requires following metadata: ``Control``, ``Assembly``, ``SemanticType``.
- ``VcfFrame[Consolidated]``
    * VcfFrame for storing target gene's consolidated variant data.
    * Requires following metadata:
- ``VcfFrame[Imported]``
    * VcfFrame for storing target gene's raw variant data.
    * Requires following metadata:
- ``VcfFrame[Phased]``
    * VcfFrame for storing target gene's phased variant data.
    * Requires following metadata:

pypgx.api.utils
===============

.. automodule:: pypgx.api.utils
   :members:
