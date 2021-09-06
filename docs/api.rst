API
***

This section describes application programming interface (API) for PyPGx.

Result and semantic type
========================

+-----+----------------------------+--------------------------------------------------------------------------------------------------------------------+
| No. | Name                       | Description                                                                                                        |
+=====+============================+====================================================================================================================+
| 1   | ``CovFrame[CopyNumber]``   | CovFrame for storing target gene's per-base copy number which is computed from read depth with control statistics. |
+-----+----------------------------+--------------------------------------------------------------------------------------------------------------------+
| 2   | ``CovFrame[ReadDepth]``    | CovFrame for storing target gene's per-base read depth which is computed from BAM files.                           |
+-----+----------------------------+--------------------------------------------------------------------------------------------------------------------+
| 3   | ``TSV[CNVCalls]``          | TSV table for storing target gene's CNV call for each sample.                                                      |
+-----+----------------------------+--------------------------------------------------------------------------------------------------------------------+
| 4   | ``TSV[Statistcs]``         | TSV table for storing control gene's various statistics on read depth.                                             |
+-----+----------------------------+--------------------------------------------------------------------------------------------------------------------+
| 5   | ``VcfFrame[Imported]``     | VcfFrame for storing target gene's raw variant data.                                                               |
+-----+----------------------------+--------------------------------------------------------------------------------------------------------------------+
| 6   | ``VcfFrame[Phased]``       | VcfFrame for storing target gene's phased variant data.                                                            |
+-----+----------------------------+--------------------------------------------------------------------------------------------------------------------+
| 7   | ``VcfFrame[Consolidated]`` | VcfFrame for storing target gene's consolidated variant data.                                                      |
+-----+----------------------------+--------------------------------------------------------------------------------------------------------------------+

pypgx.api.utils
===============

.. automodule:: pypgx.api.utils
   :members:
