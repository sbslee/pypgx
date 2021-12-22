Genes
*****

Introduction
============

This page describes gene-specific information. PyPGx currently supports a
total of 58 pharmacogenes.

Many of the genes are known to have :ref:`structural variation (SV)
<glossary:Structural variation (SV)>` including gene deletions, duplications,
and hybrids.

Some genes have a genotype-phenotype table available from the
:ref:`glossary:Clinical Pharmacogenetics Implementation Consortium (CPIC)` or
the :ref:`glossary:Pharmacogenomics Knowledge Base (PharmGKB)`. Please read
the :ref:`readme:Phenotype prediction` page for more details.

Below is a summary table:

.. list-table::
   :header-rows: 1
   :widths: 15 10 15 15 15 60

   * - Gene
     - SV
     - Phenotype
     - PharmVar
     - CPIC
     - Note
   * - ABCB1
     -
     -
     -
     -
     -
   * - :ref:`genes:CACNA1S`
     -
     - ✅
     -
     - ✅
     -
   * - :ref:`genes:CFTR`
     -
     - ✅
     -
     - ✅
     -
   * - :ref:`genes:CYP1A1`
     -
     -
     - ✅
     -
     -
   * - :ref:`genes:CYP1A2`
     -
     -
     - ✅
     -
     -
   * - :ref:`genes:CYP1B1`
     -
     -
     - ✅
     -
     -
   * - :ref:`genes:CYP2A6`
     - ✅
     -
     - ✅
     -
     - Has pseudogene (CYP2A7).
   * - :ref:`genes:CYP2A13`
     -
     -
     - ✅
     -
     -
   * - :ref:`genes:CYP2B6`
     - ✅
     - ✅
     - ✅
     - ✅
     - Has pseudogene (CYP2B7).
   * - :ref:`genes:CYP2C8`
     -
     -
     - ✅
     -
     -
   * - :ref:`genes:CYP2C9`
     -
     - ✅
     - ✅
     - ✅
     -
   * - :ref:`genes:CYP2C19`
     -
     - ✅
     - ✅
     - ✅
     -
   * - :ref:`genes:CYP2D6`
     - ✅
     - ✅
     - ✅
     - ✅
     - Has pseudogene (CYP2D7).
   * - :ref:`genes:CYP2E1`
     - ✅
     -
     - ✅
     -
     -
   * - :ref:`genes:CYP2F1`
     -
     -
     - ✅
     -
     -
   * - :ref:`genes:CYP2J2`
     -
     -
     - ✅
     -
     -
   * - :ref:`genes:CYP2R1`
     -
     -
     - ✅
     -
     -
   * - :ref:`genes:CYP2S1`
     -
     -
     - ✅
     -
     -
   * - :ref:`genes:CYP2W1`
     -
     -
     - ✅
     -
     -
   * - :ref:`genes:CYP3A4`
     -
     -
     - ✅
     -
     -
   * - :ref:`genes:CYP3A5`
     -
     - ✅
     - ✅
     - ✅
     -
   * - :ref:`genes:CYP3A7`
     -
     -
     - ✅
     -
     -
   * - :ref:`genes:CYP3A43`
     -
     -
     - ✅
     -
     -
   * - :ref:`genes:CYP4A11`
     -
     -
     - ✅
     -
     -
   * - :ref:`genes:CYP4A22`
     -
     -
     - ✅
     -
     -
   * - :ref:`genes:CYP4B1`
     -
     -
     - ✅
     -
     -
   * - :ref:`genes:CYP4F2`
     - ✅
     -
     - ✅
     -
     -
   * - :ref:`genes:CYP17A1`
     -
     -
     - ✅
     -
     -
   * - :ref:`genes:CYP19A1`
     -
     -
     - ✅
     -
     -
   * - :ref:`genes:CYP26A1`
     -
     -
     - ✅
     -
     -
   * - :ref:`genes:DPYD`
     -
     - ✅
     - ✅
     - ✅
     -
   * - :ref:`genes:F5`
     -
     - ✅
     -
     -
     -
   * - G6PD
     -
     -
     -
     -
     -
   * - :ref:`genes:GSTM1`
     - ✅
     -
     -
     -
     -
   * - GSTP1
     -
     -
     -
     -
     -
   * - :ref:`genes:GSTT1`
     - ✅
     -
     -
     -
     - Contig differs between GRCh37 and GRCh38.
   * - :ref:`genes:IFNL3`
     -
     - ✅
     -
     -
     -
   * - NAT1
     -
     -
     -
     -
     -
   * - NAT2
     -
     -
     -
     -
     -
   * - :ref:`genes:NUDT15`
     -
     - ✅
     - ✅
     - ✅
     -
   * - :ref:`genes:POR`
     -
     -
     - ✅
     -
     -
   * - :ref:`genes:PTGIS`
     -
     -
     - ✅
     -
     -
   * - :ref:`genes:RYR1`
     -
     - ✅
     - ✅
     -
     -
   * - SLC15A2
     -
     -
     -
     -
     -
   * - :ref:`genes:SLC22A2`
     - ✅
     -
     -
     -
     -
   * - :ref:`genes:SLCO1B1`
     -
     - ✅
     - ✅
     - ✅
     -
   * - SLCO1B3
     -
     -
     -
     -
     -
   * - SLCO2B1
     -
     -
     -
     -
     -
   * - :ref:`genes:SULT1A1`
     - ✅
     -
     -
     -
     -
   * - :ref:`genes:TBXAS1`
     -
     -
     - ✅
     -
     -
   * - :ref:`genes:TPMT`
     -
     - ✅
     -
     - ✅
     -
   * - :ref:`genes:UGT1A1`
     -
     - ✅
     -
     - ✅
     -
   * - :ref:`genes:UGT1A4`
     - ✅
     -
     -
     -
     -
   * - UGT2B7
     -
     -
     -
     -
     -
   * - :ref:`genes:UGT2B15`
     - ✅
     -
     -
     -
     -
   * - :ref:`genes:UGT2B17`
     - ✅
     -
     -
     -
     -
   * - VKORC1
     -
     -
     -
     - ✅
     -
   * - XPC
     -
     -
     -
     -
     -

CACNA1S
=======

Phenotype summary for CACNA1S
-----------------------------

Diplotype-phenotype mapping is used for phenotype prediction.

 .. list-table::
    :header-rows: 1

    * - Phenotype
      - Example
    * - Uncertain Susceptibility
      - Reference/Reference
    * - Malignant Hyperthermia Susceptibility
      - Reference/c.520C>T

Resources for CACNA1S
---------------------

- `Annotation of CPIC Guideline for desflurane and CACNA1S, RYR1 <https://www.pharmgkb.org/chemical/PA164749136/guidelineAnnotation/PA166180457>`__
- `CPIC® Guideline for Potent Volatile Anesthetic Agents and Succinylcholine and RYR1 and CACNA1S <https://cpicpgx.org/guidelines/cpic-guideline-for-ryr1-and-cacna1s/>`__

CFTR
====

Phenotype summary for CFTR
--------------------------

Diplotype-phenotype mapping is used for phenotype prediction.

 .. list-table::
    :header-rows: 1

    * - Phenotype
      - Example
    * - Favorable Response
      - Reference/G551D
    * - Unfavorable Response
      - F508del/F508del
    * - Indeterminate
      - Reference/F508del

Resources for CFTR
------------------

- `Annotation of CPIC Guideline for ivacaftor and CFTR <https://www.pharmgkb.org/chemical/PA165950341/guidelineAnnotation/PA166114461>`__
- `CPIC® Guideline for Ivacaftor and CFTR <https://cpicpgx.org/guidelines/guideline-for-ivacaftor-and-cftr/>`__

CYP1A1
======

Resources for CYP1A1
--------------------

- `PharmVar CYP1A1 page <https://www.pharmvar.org/gene/CYP1A1>`__

CYP1A2
======

Resources for CYP1A2
--------------------

- `PharmVar CYP1A2 page <https://www.pharmvar.org/gene/CYP1A2>`__

CYP1B1
======

Resources for CYP1B1
--------------------

- `PharmVar CYP1B1 page <https://www.pharmvar.org/gene/CYP1B1>`__

CYP2A6
======

SV summary for CYP2A6
---------------------

Below is comprehensive summary of SV described from real NGS studies:

.. list-table::
   :header-rows: 1

   * - Star Allele
     - SV Name
     - Genotype
     - Reference
     - GRCh37
     - GRCh38
     - Data Type
     - Source
     - Coriell ID
   * - \*4
     - DeletionHet
     - \*1/\*4
     -
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-1.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-1.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA18617
   * - \*4
     - DeletionHom
     - \*4/\*4
     -
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-2.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-2.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA18952
   * - \*1x2
     - Duplication
     - \*1x2/\*25
     - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-3.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-3.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA18861
   * - Indeterminate
     - Hybrid
     - Indeterminate
     -
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-4.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-4.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - HG00436

Filtered alleles for CYP2A6
---------------------------

Some alleles in PharmVar will not be called by PyPGx because one or more of their variants have a high false positive rate, likely due to read misalignment to the *CYP2A7* pseudogene. Those alleles are listed in below table. If problematic variants are present in gnomAD, their links are provided so that you can look at filtering status, allele imbalance for heterozygotes, etc.

.. list-table::
   :widths: 25 25 25 25
   :header-rows: 1

   * - Problematic Variant
     - Star Alleles
     - GRCh37
     - GRCh38
   * - rs143731390 (N438Y)
     - \*35
     - `22-42523514-C-T <https://gnomad.broadinstitute.org/variant/19-41349874-T-A?dataset=gnomad_r2_1>`__
     - `22-42127512-C-T <https://gnomad.broadinstitute.org/variant/19-40843969-T-A?dataset=gnomad_r3>`__

CYP2A13
=======

Resources for CYP2A13
---------------------

- `PharmVar CYP2A13 page <https://www.pharmvar.org/gene/CYP2A13>`__

CYP2B6
======

SV summary for CYP2B6
---------------------

Below is comprehensive summary of SV described from real NGS studies:

.. list-table::
   :header-rows: 1

   * - Star Allele
     - SV Name
     - Genotype
     - Reference
     - GRCh37
     - GRCh38
     - Data Type
     - Source
     - Coriell ID
   * - \*29
     - Hybrid
     - \*6/\*29
     -
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2B6-1.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2B6-1.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA19178

Phenotype summary for CYP2B6
----------------------------

Diplotype-phenotype mapping is used for phenotype prediction.

 .. list-table::
    :header-rows: 1

    * - Phenotype
      - Example
    * - Ultrarapid Metabolizer
      - \*4/\*4
    * - Rapid Metabolizer
      - \*1/\*4
    * - Normal Metabolizer
      - \*1/\*2
    * - Intermediate Metabolizer
      - \*1/\*29
    * - Poor Metabolizer
      - \*6/\*6
    * - Indeterminate
      - \*1/\*3

Resources for CYP2B6
--------------------

- `CPIC® Guideline for Efavirenz based on CYP2B6 genotype <https://cpicpgx.org/guidelines/cpic-guideline-for-efavirenz-based-on-cyp2b6-genotype/>`__

CYP2C8
======

Resources for CYP2C8
--------------------

- `PharmVar CYP2C8 page <https://www.pharmvar.org/gene/CYP2C8>`__

CYP2C9
======

Phenotype summary for CYP2C9
----------------------------

Activity score is used for phenotype prediction.

.. list-table::
   :header-rows: 1

   * - Phenotype
     - Activity Score
     - Example
   * - Normal Metabolizer
     - 2 == score
     - \*1/\*1
   * - Intermediate Metabolizer
     - 1 <= score < 2
     - \*1/\*2
   * - Poor Metabolizer
     - 0 <= score < 1
     - \*2/\*3

Resources for CYP2C9
--------------------

- `CPIC® Guideline for NSAIDs based on CYP2C9 genotype <https://cpicpgx.org/guidelines/cpic-guideline-for-nsaids-based-on-cyp2c9-genotype/>`__

CYP2C19
=======

Phenotype summary for CYP2C19
-----------------------------

Diplotype-phenotype mapping is used for phenotype prediction.

.. list-table::
   :header-rows: 1

   * - Phenotype
     - Example
   * - Ultrarapid Metabolizer
     - \*17/\*17
   * - Rapid Metabolizer
     - \*1/\*17
   * - Normal Metabolizer
     - \*1/\*1
   * - Likely Intermediate Metabolizer
     - \*1/\*10
   * - Intermediate Metabolizer
     - \*1/\*2
   * - Likely Poor Metabolizer
     - \*10/\*22
   * - Poor Metabolizer
     - \*2/\*2
   * - Indeterminate
     - \*1/\*12

Resources for CYP2C19
---------------------

- `CPIC® Guideline for Voriconazole and CYP2C19 <https://cpicpgx.org/guidelines/guideline-for-voriconazole-and-cyp2c19/>`__

CYP2D6
======

SV summary for CYP2D6
---------------------

Below is comprehensive summary of SV described from real NGS studies:

.. list-table::
   :header-rows: 1

   * - Star Allele
     - SV Name
     - Genotype
     - Reference
     - GRCh37
     - GRCh38
     - Data Type
     - Source
     - Coriell ID
   * - \*5
     - DeletionHet
     - \*5/\*29
     -
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-1.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-1.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA18861
   * - \*5
     - DeletionHom
     - \*5/\*5
     -
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-6.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-6.png>`
     - WGS
     -
     -
   * - \*4x2
     - Duplication
     - \*2/\*4x2
     -
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-2.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-2.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA19819
   * - \*68+\*4
     - Tandem1
     - \*139/\*68+\*4
     -
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-3.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-3.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA11832
   * - \*36+\*10
     - Tandem2A
     - \*2/\*36+\*10
     -
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-4.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-4.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA18564
   * - \*36x2+\*10
     - Tandem2B
     - \*1/\*36x2+\*10
     -
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-5.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-5.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA18524
   * - \*36x3+\*10
     - Tandem2C
     - \*1/\*36x3+\*10
     -
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-7.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-7.png>`
     - WGS
     -
     -

Phenotype summary for CYP2D6
----------------------------

Activity score is used for phenotype prediction.

.. list-table::
   :header-rows: 1

   * - Phenotype
     - Activity Score
     - Example
   * - Ultrarapid Metabolizer
     - 2.5 <= score
     - \*1/\*2x2
   * - Normal Metabolizer
     - 1.25 <= score < 2.5
     - \*1/\*1
   * - Intermediate Metabolizer
     - 0.25 <= score < 1.25
     - \*1/\*4
   * - Poor Metabolizer
     - 0 <= score < 0.25
     - \*4/\*5

Filtered alleles for CYP2D6
---------------------------

Some alleles in PharmVar will not be called by PyPGx because one or more of their variants have a high false positive rate, likely due to read misalignment to the *CYP2D7* pseudogene. Those alleles are listed in below table. If problematic variants are present in gnomAD, their links are provided so that you can look at filtering status, allele imbalance for heterozygotes, etc.

.. list-table::
   :widths: 25 25 25 25
   :header-rows: 1

   * - Problematic Variant
     - Star Alleles
     - GRCh37
     - GRCh38
   * - rs769157652 (E410K)
     - \*27, \*32
     - `22-42522940-C-T <https://gnomad.broadinstitute.org/variant/22-42522940-C-T?dataset=gnomad_r2_1>`__
     - `22-42126938-C-T <https://gnomad.broadinstitute.org/variant/22-42126938-C-T?dataset=gnomad_r3>`__
   * - rs61745683 (V370I)
     - \*122
     - `22-42523514-C-T <https://gnomad.broadinstitute.org/variant/22-42523514-C-T?dataset=gnomad_r2_1>`__
     - `22-42127512-C-T <https://gnomad.broadinstitute.org/variant/22-42127512-C-T?dataset=gnomad_r3>`__
   * - rs1058172 (R365H)
     - \*139
     - `22-42523528-C-T <https://gnomad.broadinstitute.org/variant/22-42523528-C-T?dataset=gnomad_r2_1>`__
     - `22-42127526-C-T <https://gnomad.broadinstitute.org/variant/22-42127526-C-T?dataset=gnomad_r3>`__
   * - rs202102799 (Y355C)
     - \*127
     - `22-42523558-T-C <https://gnomad.broadinstitute.org/variant/22-42523558-T-C?dataset=gnomad_r2_1>`__
     - `22-42127556-T-C <https://gnomad.broadinstitute.org/variant/22-42127556-T-C?dataset=gnomad_r3>`__
   * - rs17002853 (L231P)
     - \*131
     - `22-42524327-A-G <https://gnomad.broadinstitute.org/variant/22-42524327-A-G?dataset=gnomad_r2_1>`__
     - `22-42128325-A-G <https://gnomad.broadinstitute.org/variant/22-42128325-A-G?dataset=gnomad_r3>`__

Resources for CYP2D6
--------------------

- `CPIC® Guideline for Tamoxifen based on CYP2D6 genotype <https://cpicpgx.org/guidelines/cpic-guideline-for-tamoxifen-based-on-cyp2d6-genotype/>`__

CYP2E1
======

SV summary for CYP2E1
---------------------

Below is comprehensive summary of SV described from real NGS studies:

.. list-table::
   :header-rows: 1

   * - Star Allele
     - SV Name
     - Genotype
     - Reference
     - GRCh37
     - GRCh38
     - Data Type
     - Source
     - Coriell ID
   * - \*S1
     - PartialDuplication
     - \*1/\*S1
     - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2E1-1.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2E1-1.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA19920
   * - \*1x2
     - Duplication
     - \*1/\*1x2
     -
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2E1-4.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2E1-4.png>`
     - WGS
     -
     -
   * - \*7x2
     - Duplication
     - \*1/\*7x2
     - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2E1-2.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2E1-2.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA19095
   * - \*7x3
     - Multiplication
     - \*7/\*7x3
     -
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2E1-3.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2E1-3.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA19908

CYP2F1
======

Resources for CYP2F1
--------------------

- `PharmVar CYP2F1 page <https://www.pharmvar.org/gene/CYP2F1>`__

CYP2J2
======

Resources for CYP2J2
--------------------

- `PharmVar CYP2J2 page <https://www.pharmvar.org/gene/CYP2J2>`__

CYP2R1
======

Resources for CYP2R1
--------------------

- `PharmVar CYP2R1 page <https://www.pharmvar.org/gene/CYP2R1>`__

CYP2S1
======

Resources for CYP2S1
--------------------

- `PharmVar CYP2S1 page <https://www.pharmvar.org/gene/CYP2S1>`__

CYP2W1
======

Resources for CYP2W1
--------------------

- `PharmVar CYP2W1 page <https://www.pharmvar.org/gene/CYP2W1>`__

CYP3A4
======

Resources for CYP3A4
--------------------

- `PharmVar CYP3A4 page <https://www.pharmvar.org/gene/CYP3A4>`__

CYP3A5
======

Phenotype summary for CYP3A5
----------------------------

Diplotype-phenotype mapping is used for phenotype prediction.

.. list-table::
   :header-rows: 1

   * - Phenotype
     - Example
   * - Normal Metabolizer
     - \*1/\*1
   * - Intermediate Metabolizer
     - \*1/\*3
   * - Possible Intermediate Metabolizer
     - \*1/\*2
   * - Poor Metabolizer
     - \*6/\*6
   * - Indeterminate
     - \*2/\*2

Resources for CYP3A5
--------------------

- `CPIC® Guideline for Tacrolimus and CYP3A5 <https://cpicpgx.org/guidelines/guideline-for-tacrolimus-and-cyp3a5/>`__

CYP3A7
======

Resources for CYP3A7
--------------------

- `PharmVar CYP3A7 page <https://www.pharmvar.org/gene/CYP3A7>`__

CYP3A43
=======

Resources for CYP3A43
---------------------

- `PharmVar CYP3A43 page <https://www.pharmvar.org/gene/CYP3A43>`__

CYP4A11
=======

Resources for CYP4A11
---------------------

- `PharmVar CYP4A11 page <https://www.pharmvar.org/gene/CYP4A11>`__

CYP4A22
=======

Resources for CYP4A22
---------------------

- `PharmVar CYP4A22 page <https://www.pharmvar.org/gene/CYP4A22>`__

CYP4B1
======

Resources for CYP4B1
--------------------

- `PharmVar CYP4B1 page <https://www.pharmvar.org/gene/CYP4B1>`__

CYP4F2
======

SV summary for CYP4F2
---------------------

Below is comprehensive summary of SV described from real NGS studies:

.. list-table::
  :header-rows: 1

  * - Star Allele
    - SV Name
    - Genotype
    - Reference
    - GRCh37
    - GRCh38
    - Data Type
    - Source
    - Coriell ID
  * - \*DEL
    - DeletionHet
    - \*1/\*DEL
    -
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP4F2-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP4F2-1.png>`
    - WGS
    - 
    -

Resources for CYP4F2
--------------------

- `PharmVar CYP4F2 page <https://www.pharmvar.org/gene/CYP4F2>`__

CYP17A1
=======

Resources for CYP17A1
---------------------

- `PharmVar CYP17A1 page <https://www.pharmvar.org/gene/CYP17A1>`__

CYP19A1
=======

Resources for CYP19A1
---------------------

- `PharmVar CYP19A1 page <https://www.pharmvar.org/gene/CYP19A1>`__

CYP26A1
=======

Resources for CYP26A1
---------------------

- `PharmVar CYP26A1 page <https://www.pharmvar.org/gene/CYP26A1>`__

DPYD
====

Phenotype summary for DPYD
--------------------------

Activity score is used for phenotype prediction.

.. list-table::
   :header-rows: 1

   * - Phenotype
     - Activity Score
     - Example
   * - Normal Metabolizer
     - 2 == score
     - Reference/Reference
   * - Intermediate Metabolizer
     - 1 <= score < 2
     - Reference/c.1905+1G>A (\*2A)
   * - Poor Metabolizer
     - 0 <= score < 1
     - c.295_298delTCAT (\*7)/c.703C>T (\*8)

Resources for DPYD
------------------

- `CPIC® Guideline for Fluoropyrimidines and DPYD <https://cpicpgx.org/guidelines/guideline-for-fluoropyrimidines-and-dpyd/>`__

F5
==

Phenotype summary for F5
------------------------

Diplotype-phenotype mapping is used for phenotype prediction.

 .. list-table::
    :header-rows: 1

    * - Phenotype
      - Example
    * - Favorable Response
      - Reference/Reference
    * - Unfavorable Response
      - Reference/Leiden

Resources for F5
----------------

- `Annotation of DPWG Guideline for hormonal contraceptives for systemic use and F5 <https://www.pharmgkb.org/chemical/PA452637/guidelineAnnotation/PA166104955>`__

GSTM1
=====

SV summary for GSTM1
--------------------

This gene is known to have an extremely high rate of gene deletion
polymorphism in the population and thus requires SV analysis.

Below is comprehensive summary of SV described from real NGS studies:

.. list-table::
  :header-rows: 1

  * - Star Allele
    - SV Name
    - Genotype
    - Reference
    - GRCh37
    - GRCh38
    - Data Type
    - Source
    - Coriell ID
  * - \*0
    - DeletionHet
    - \*0/\*A
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTM1-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTM1-1.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA18855
  * - \*0
    - DeletionHom
    - \*0/\*0
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTM1-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTM1-2.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA10831
  * - \*Ax2
    - Duplication
    - \*A/\*Ax2
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTM1-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTM1-3.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA19908
  * - \*Bx2
    - Duplication
    - \*A/\*Bx2
    -
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTM1-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTM1-4.png>`
    - WGS
    -
    -

GSTT1
=====

GRCh38 data for GSTT1
---------------------

*GSTT1* is located on ``chr22`` for GRCh37 but on ``chr22_KI270879v1_alt``
for GRCh38. Therefore, if you are interested in genotyping this gene with
GRCh38 data, then you must have sequence reads mapped to the ALT contig.
For more details, please read the :ref:`readme:GRCh37 vs. GRCh38` page.

SV summary for GSTT1
--------------------

This gene is known to have an extremely high rate of gene deletion
polymorphism in the population and thus requires SV analysis.

Below is comprehensive summary of SV described from real NGS studies:

.. list-table::
  :header-rows: 1

  * - Star Allele
    - SV Name
    - Genotype
    - Reference
    - GRCh37
    - GRCh38
    - Data Type
    - Source
    - Coriell ID
  * - \*0
    - DeletionHet
    - \*0/\*A
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTT1-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTT1-1.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA19908
  * - \*0
    - DeletionHom
    - \*0/\*0
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTT1-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTT1-2.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA11832

IFNL3
=====

Phenotype summary for IFNL3
---------------------------

Diplotype-phenotype mapping is used for phenotype prediction.

.. list-table::
   :header-rows: 1

   * - Phenotype
     - Example
   * - Favorable Response
     - Reference/Reference
   * - Unfavorable Response
     - Reference/rs12979860
   * - Indeterminate
     - Reference/rs8099917

Resources for IFNL3
-------------------

- `Annotation of CPIC Guideline for peginterferon alfa-2a,peginterferon alfa-2b,ribavirin and IFNL3 <https://www.pharmgkb.org/guidelineAnnotation/PA166110235>`__

NUDT15
======

Phenotype summary for NUDT15
----------------------------

Diplotype-phenotype mapping is used for phenotype prediction.

.. list-table::
   :header-rows: 1

   * - Phenotype
     - Example
   * - Normal Metabolizer
     - \*1/\*1
   * - Intermediate Metabolizer
     - \*1/\*2
   * - Possible Intermediate Metabolizer
     - \*3/\*4
   * - Poor Metabolizer
     - \*2/\*3
   * - Indeterminate
     - \*1/\*4

Resources for NUDT15
--------------------

- `CPIC® Guideline for Thiopurines and TPMT and NUDT15 <https://cpicpgx.org/guidelines/guideline-for-thiopurines-and-tpmt/>`__

POR
===

Resources for POR
-----------------

- `PharmVar POR page <https://www.pharmvar.org/gene/POR>`__

PTGIS
=====

Resources for PTGIS
-------------------

- `PharmVar PTGIS page <https://www.pharmvar.org/gene/PTGIS>`__

RYR1
====

Phenotype summary for RYR1
--------------------------

Diplotype-phenotype mapping is used for phenotype prediction.

.. list-table::
  :header-rows: 1

  * - Phenotype
    - Example
  * - Uncertain Susceptibility
    - Reference/Reference
  * - Malignant Hyperthermia Susceptibility
    - Reference/c.103T>C

Resources for RYR1
------------------

- `Annotation of CPIC Guideline for desflurane and CACNA1S, RYR1 <https://www.pharmgkb.org/chemical/PA164749136/guidelineAnnotation/PA166180457>`__
- `CPIC® Guideline for Potent Volatile Anesthetic Agents and Succinylcholine and RYR1 and CACNA1S <https://cpicpgx.org/guidelines/cpic-guideline-for-ryr1-and-cacna1s/>`__

SLC22A2
=======

SV summary for SLC22A2
----------------------

Below is comprehensive summary of SV described from real NGS studies:

.. list-table::
  :header-rows: 1

  * - Star Allele
    - SV Name
    - Genotype
    - Reference
    - GRCh37
    - GRCh38
    - Data Type
    - Source
    - Coriell ID
  * - \*S1
    - Intron9Deletion
    - \*1/\*S1
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SLC22A2-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SLC22A2-1.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA18855
  * - \*S2
    - Exon11Deletion
    - \*1/\*S2
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SLC22A2-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SLC22A2-2.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA19819

SLCO1B1
=======

Phenotype summary for SLCO1B1
-----------------------------

Diplotype-phenotype mapping is used for phenotype prediction.

.. list-table::
   :header-rows: 1

   * - Phenotype
     - Example
   * - Possible Increased Function
     - \*1A/\*35
   * - Normal Function
     - \*1A/\*1A
   * - Possible Decreased Function
     - \*1A/\*2
   * - Decreased Function
     - \*1A/\*5
   * - Possible Poor Function
     - \*2/\*2
   * - Poor Function
     - \*5/\*5
   * - Indeterminate
     - \*1A/\*7

Resources for SLCO1B1
---------------------

- `CPIC® Guideline for Simvastatin and SLCO1B1 <https://cpicpgx.org/guidelines/guideline-for-simvastatin-and-slco1b1/>`__

SULT1A1
=======

SV summary for SULT1A1
----------------------

Below is comprehensive summary of SV described from real NGS studies:

.. list-table::
  :header-rows: 1

  * - Star Allele
    - SV Name
    - Genotype
    - Reference
    - GRCh37
    - GRCh38
    - Data Type
    - Source
    - Coriell ID
  * - \*DEL
    - DeletionHet
    - \*1/\*DEL
    -
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SULT1A1-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SULT1A1-1.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA18942
  * - \*1x2
    - Duplication
    - \*1x2/\*2
    -
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SULT1A1-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SULT1A1-2.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA18509
  * - \*1x3
    - Multiplication1
    - \*1x3/\*2
    -
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SULT1A1-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SULT1A1-3.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA18868
  * - \*1x4
    - Multiplication2
    - \*1x4/\*2
    -
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SULT1A1-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SULT1A1-4.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA18484
  * - \*1x3, \*2x2
    - Multiplication2
    - \*1x3/\*2x2
    -
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SULT1A1-5.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SULT1A1-5.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA19143

TBXAS1
======

Resources for TBXAS1
--------------------

- `PharmVar TBXAS1 page <https://www.pharmvar.org/gene/TBXAS1>`__

TPMT
====

Phenotype summary for TPMT
--------------------------

Diplotype-phenotype mapping is used for phenotype prediction.

.. list-table::
   :header-rows: 1

   * - Phenotype
     - Example
   * - Normal Metabolizer
     - \*1/\*1
   * - Possible Intermediate Metabolizer
     - \*3A/\*12
   * - Intermediate Metabolizer
     - \*1/\*2
   * - Poor Metabolizer
     - \*2/\*3A
   * - Indeterminate
     - \*1/\*18

Resources for TPMT
------------------

- `CPIC® Guideline for Thiopurines and TPMT and NUDT15 <https://cpicpgx.org/guidelines/guideline-for-thiopurines-and-tpmt/>`__

UGT1A1
======

Phenotype summary for UGT1A1
----------------------------

Diplotype-phenotype mapping is used for phenotype prediction.

.. list-table::
   :header-rows: 1

   * - Phenotype
     - Example
   * - Normal Metabolizer
     - \*1/\*1
   * - Intermediate Metabolizer
     - \*1/\*6
   * - Poor Metabolizer
     - \*6/\*27
   * - Indeterminate
     - \*28/\*80

Resources for UGT1A1
--------------------

- `CPIC® Guideline for Atazanavir and UGT1A1 <https://cpicpgx.org/guidelines/guideline-for-atazanavir-and-ugt1a1/>`__

UGT1A4
======

SV summary for UGT1A4
---------------------

Below is comprehensive summary of SV described from real NGS studies:

.. list-table::
  :header-rows: 1

  * - Star Allele
    - SV Name
    - Genotype
    - Reference
    - GRCh37
    - GRCh38
    - Data Type
    - Source
    - Coriell ID
  * - \*S1
    - Intron1DeletionA
    - \*1/\*S1
    -
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT1A4-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT1A4-1.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA19908
  * - \*S2
    - Intron1DeletionB
    - \*1/\*S2
    -
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT1A4-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT1A4-2.png>`
    - WGS
    -
    -

UGT2B15
=======

SV summary for UGT2B15
----------------------

Below is comprehensive summary of SV described from real NGS studies:

.. list-table::
  :header-rows: 1

  * - Star Allele
    - SV Name
    - Genotype
    - Reference
    - GRCh37
    - GRCh38
    - Data Type
    - Source
    - Coriell ID
  * - \*S1
    - PartialDeletion
    - \*4/\*S1
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT2B15-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT2B15-1.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA11993

UGT2B17
=======

SV summary for UGT2B17
----------------------

This gene is known to have an extremely high rate of gene deletion
polymorphism in the population and thus requires SV analysis.

Below is comprehensive summary of SV described from real NGS studies:

.. list-table::
  :header-rows: 1

  * - Star Allele
    - SV Name
    - Genotype
    - Reference
    - GRCh37
    - GRCh38
    - Data Type
    - Source
    - Coriell ID
  * - \*2
    - DeletionHet
    - \*1/\*2
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT2B17-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT2B17-1.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA18855
  * - \*2
    - DeletionHom
    - \*2/\*2
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT2B17-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT2B17-2.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA18617
