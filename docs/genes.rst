Genes
*****

This page describes gene-specific information. PyPGx currently supports
genotyping of a total of 88 pharmacogenes.

In order to provide the most accurate information, this page borrows heavily
from the works of the :ref:`glossary:Clinical Pharmacogenetics Implementation
Consortium (CPIC)` and the :ref:`glossary:Pharmacogenomics Knowledge Base
(PharmGKB)`. All curated contents published by CPIC and PharmGKB are
available free of restriction under the `CC0 1.0 Universal (CC0 1.0) Public
Domain Dedication <https://cpicpgx.org/license/>`__ and the `Creative Commons
Attribution-ShareAlike 4.0 International License <https://www.pharmgkb.org/
page/dataUsagePolicy>`__, repsectively.

Many of the genes are known to have :ref:`structural variation (SV)
<glossary:Structural variation (SV)>` including gene deletions, duplications,
and hybrids. Please read the :ref:`readme:Structural variation detection`
page for more details.

Some genes have a genotype-phenotype table available from CPIC or PharmGKB.
Please read the :ref:`readme:Phenotype prediction` page for more details.

Below is a summary table:

.. list-table::
   :header-rows: 1

   * - Gene
     - Variants
     - SV
     - Phenotype
     - PharmVar
     - CPIC
     - Function
     - GRCh37
     - GRCh38
     - Notes
   * - ABCB1
     - ✅
     -
     -
     -
     -
     - Disposition
     - `chr7:87130178-87345639 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A87130178%2D87345639&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr7:87500862-87716323 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A87500862%2D87716323&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:ABCG2`
     - ✅
     -
     - ✅
     -
     - ✅
     - Disposition
     - `chr4:89008420-89082791 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr4%3A89008420%2D89082791&hgsid=1298429733_BSyanNFtoxsgwNZmMlPdvfYamJmW>`__
     - `chr4:88087268-88161639 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr4%3A88087268%2D88161639&hgsid=1298429733_BSyanNFtoxsgwNZmMlPdvfYamJmW>`__
     -
   * - ACYP2
     - ✅
     -
     -
     -
     -
     - Metabolism
     - `chr2:54339751-54535437 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2%3A54339751%2D54535437&hgsid=2075643426_jPVGUWohm3uLIv8bdQ6vCs7HSa7v>`__
     - `chr2:54112614-54308300 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2%3A54112614%2D54308300&hgsid=2075643426_jPVGUWohm3uLIv8bdQ6vCs7HSa7v>`__
     - NM_138448.4 was used as the main transcript.
   * - ADRA2A
     - ✅
     -
     -
     -
     -
     - Absorption
     - `chr10:112833786-112843665 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A112833786%2D112843665&hgsid=2091670486_9UK3D1ZLOfc1KFoEhawm0PTpFAer>`__
     - `chr10:111074028-111083907 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A111074028%2D111083907&hgsid=2091670486_9UK3D1ZLOfc1KFoEhawm0PTpFAer>`__
     - NM_000681.4 was used as the main transcript.
   * - ADRB2
     - ✅
     -
     -
     -
     -
     - Target
     - `chr5:148203173-148211186 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr5%3A148203173%2D148211186&hgsid=2092325066_K6VvMJMxdGa5dur0C6kdxfWOxiwZ>`__
     - `chr5:148823610-148831623 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr5%3A148823610%2D148831623&hgsid=2092325066_K6VvMJMxdGa5dur0C6kdxfWOxiwZ>`__
     - NM_000024.6 was used as the main transcript.
   * - ANKK1
     - ✅
     -
     -
     -
     -
     - Metabolism
     - `chr11:113255500-113274138 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr11%3A113255500%2D113274138&hgsid=2097319032_Yh0JY7nKAHHYsaU4mObJ1Km7fTXg>`__
     - `chr11:113384778-113403416 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr11%3A113384778%2D113403416&hgsid=2097319032_Yh0JY7nKAHHYsaU4mObJ1Km7fTXg>`__
     - NM_178510.2 was used as the main transcript.
   * - APOE
     - ✅
     -
     -
     -
     -
     - Efficacy
     - `chr19:45406052-45415650 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A45406052%2D45415650&hgsid=2101995320_cUhlNCjppaGJN0dqC3j4Rz9XAqaD>`__
     - `chr19:44902795-44912393 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A44902795%2D44912393&hgsid=2101995320_cUhlNCjppaGJN0dqC3j4Rz9XAqaD>`__
     - NM_000041.4 was used as the main transcript.
   * - ATM
     - ✅
     -
     -
     -
     -
     - Efficacy
     - `chr11:108090793-108286161 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr11%3A108090793%2D108286161&hgsid=2101995320_cUhlNCjppaGJN0dqC3j4Rz9XAqaD>`__
     - `chr11:108220066-108415434 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr11%3A108220066%2D108415434&hgsid=2101995320_cUhlNCjppaGJN0dqC3j4Rz9XAqaD>`__
     - NM_000051.4 was used as the main transcript.
   * - BCHE
     - ✅
     -
     -
     -
     -
     - Efficacy
     - `chr3:165487691-165558211 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr3%3A165487691%2D165558211&hgsid=2101995320_cUhlNCjppaGJN0dqC3j4Rz9XAqaD>`__
     - `chr3:165769903-165840423 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr3%3A165769903%2D165840423&hgsid=2101995320_cUhlNCjppaGJN0dqC3j4Rz9XAqaD>`__
     - NM_000055.4 was used as the main transcript.
   * - BDNF
     - ✅
     -
     -
     -
     -
     - Efficacy
     - `chr11:27673439-27725002 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr11%3A27673439%2D27725002&hgsid=2104479616_zVjvFDwop6R0m6VrmJkHXmy2Wayd>`__
     - `chr11:27651892-27703455 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr11%3A27651892%2D27703455&hgsid=2104479616_zVjvFDwop6R0m6VrmJkHXmy2Wayd>`__
     - NM_001709.5 was used as the main transcript.
   * - :ref:`genes:CACNA1S`
     - ✅
     -
     - ✅
     -
     - ✅
     - Target
     - `chr1:201005639-201084694 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A201005639%2D201084694&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr1:201036511-201115426 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A201036511%2D201115426&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CFTR`
     - ✅
     -
     - ✅
     -
     - ✅
     - Target
     - `chr7:117117016-117311719 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A117117016%2D117311719&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr7:117477024-117671665 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A117477024%2D117671665&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - COMT
     - ✅
     -
     -
     -
     -
     - Metabolism
     - `chr22:19917359-19962473 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr22%3A19917359%2D19962473&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr22:19929836-19974950 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr22%3A19929836%2D19974950&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP1A1`
     - ✅
     -
     -
     - ✅
     -
     - Metabolism
     - `chr15:75008882-75020951 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr15%3A75008882%2D75020951&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr15:74716541-74728528 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr15%3A74716541%2D74728528&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP1A2`
     - ✅
     -
     -
     - ✅
     -
     - Metabolism
     - `chr15:75038183-75051941 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr15%3A75038183%2D75051941&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr15:74745844-74759607 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr15%3A74745844%2D74759607&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP1B1`
     - ✅
     -
     -
     - ✅
     -
     - Metabolism
     - `chr2:38291745-38306323 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2%3A38291745%2D38306323&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr2:38064602-38079181 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2%3A38064602%2D38079181&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP2A6`
     - ✅
     - ✅
     -
     - ✅
     -
     - Metabolism
     - `chr19:41339442-41396352 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A41339442%2D41396352&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr19:40833540-40890447 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A40833540%2D40890447&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - CYP2A6 has pseudogene (CYP2A7).
   * - :ref:`genes:CYP2A13`
     - ✅
     -
     -
     - ✅
     -
     - Metabolism
     - `chr19:41574355-41622100 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A41574355%2D41622100&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr19:41068450-41116195 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A41068450%2D41116195&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP2B6`
     - ✅
     - ✅
     - ✅
     - ✅
     - ✅
     - Metabolism
     - `chr19:41427203-41534301 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A41427203%2D41534301&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr19:40921281-41028398 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A40921281%2D41028398&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - CYP2B6 has pseudogene (CYP2B7).
   * - :ref:`genes:CYP2C8`
     - ✅
     -
     -
     - ✅
     -
     - Metabolism
     - `chr10:96793528-96832254 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A96793528%2D96832254&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr10:95033771-95072497 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A95033771%2D95072497&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP2C9`
     - ✅
     -
     - ✅
     - ✅
     - ✅
     - Metabolism
     - `chr10:96695414-96752148 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A96695414%2D96752148&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr10:94935657-94993091 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A94935657%2D94993091&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP2C19`
     - ✅
     -
     - ✅
     - ✅
     - ✅
     - Metabolism
     - `chr10:96519437-96615962 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A96519437%2D96615962&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr10:94759680-94858547 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A94759680%2D94858547&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP2D6`
     - ✅
     - ✅
     - ✅
     - ✅
     - ✅
     - Metabolism
     - `chr22:42512500-42551883 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr22%3A42512500%2D42551883&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr22:42116498-42155810 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr22%3A42116498%2D42155810&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - CYP2D6 has pseudogene (CYP2D7).
   * - :ref:`genes:CYP2E1`
     - ✅
     - ✅
     -
     - ✅
     -
     - Metabolism
     - `chr10:135330866-135362620 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A135330866%2D135362620&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr10:133517362-133549123 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A133517362%2D133549123&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP2F1`
     - ✅
     -
     -
     - ✅
     -
     - Metabolism
     - `chr19:41617336-41637286 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A41617336%2D41637286&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr19:41111431-41131381 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A41111431%2D41131381&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP2J2`
     - ✅
     -
     -
     - ✅
     -
     - Metabolism
     - `chr1:60355979-60395470 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A60355979%2D60395470&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr1:59890307-59929773 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A59890307%2D59929773&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP2R1`
     - ✅
     -
     -
     - ✅
     -
     - Metabolism
     - `chr11:14896554-14916751 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr11%3A14896554%2D14916751&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr11:14875008-14895205 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr11%3A14875008%2D14895205&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP2S1`
     - ✅
     -
     -
     - ✅
     -
     - Metabolism
     - `chr19:41696111-41716444 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A41696111%2D41716444&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr19:41190218-41210539 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A41190218%2D41210539&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP2W1`
     - ✅
     -
     -
     - ✅
     -
     - Metabolism
     - `chr7:1019834-1032276 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A1019834%2D1032276&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr7:980180-992640 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A980180%2D992640&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP3A4`
     - ✅
     -
     -
     - ✅
     -
     - Metabolism
     - `chr7:99351582-99384811 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A99351582%2D99384811&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr7:99753966-99787184 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A99753966%2D99787184&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP3A5`
     - ✅
     -
     - ✅
     - ✅
     - ✅
     - Metabolism
     - `chr7:99242811-99280649 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A99242811%2D99280649&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr7:99645193-99682996 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A99645193%2D99682996&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP3A7`
     - ✅
     -
     -
     - ✅
     -
     - Metabolism
     - `chr7:99299659-99335823 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A99299659%2D99335823&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr7:99702035-99738196 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A99702035%2D99738196&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP3A43`
     - ✅
     -
     -
     - ✅
     -
     - Metabolism
     - `chr7:99422635-99466727 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A99422635%2D99466727&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr7:99825012-99869093 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A99825012%2D99869093&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP4A11`
     - ✅
     -
     -
     - ✅
     -
     - Metabolism
     - `chr1:47391859-47410148 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A47391859%2D47410148&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr1:46926187-46944476 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A46926187%2D46944476&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP4A22`
     - ✅
     -
     -
     - ✅
     -
     - Metabolism
     - `chr1:47600112-47618399 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A47600112%2D47618399&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr1:47134440-47152727 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A47134440%2D47152727&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP4B1`
     - ✅
     -
     -
     - ✅
     -
     - Metabolism
     - `chr1:47261669-47288021 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A47261669%2D47288021&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr1:46796045-46822413 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A46796045%2D46822413&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP4F2`
     - ✅
     - ✅
     -
     - ✅
     -
     - Metabolism
     - `chr19:15973833-16023930 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A15973833%2D16023930&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr19:15863022-15913074 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A15863022%2D15913074&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP17A1`
     - ✅
     -
     -
     - ✅
     -
     - Metabolism
     - `chr10:104587287-104600170 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A104587287%2D104600170&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr10:102827530-102840413 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A102827530%2D102840413&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP19A1`
     - ✅
     -
     -
     - ✅
     -
     - Metabolism
     - `chr15:51497253-51633795 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr15%3A51497253%2D51633795&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr15:51205056-51341596 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr15%3A51205056%2D51341596&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP26A1`
     - ✅
     -
     -
     - ✅
     -
     - Metabolism
     - `chr10:94830646-94840641 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A94830646%2D94840641&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr10:93070892-93080885 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A93070892%2D93080885&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - DBH
     - ✅
     -
     -
     -
     -
     - Efficacy
     - `chr9:136498484-136527451 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr9%3A136498484%2D136527451&hgsid=2104624522_XC5vAwAl5AD08pCV2WxGd1Z49TF6>`__
     - `chr9:133633362-133662329 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr9%3A133633362%2D133662329&hgsid=2104624522_XC5vAwAl5AD08pCV2WxGd1Z49TF6>`__
     - NM_000787.4 was used as the main transcript.
   * - :ref:`genes:DPYD`
     - ✅
     -
     - ✅
     - ✅
     - ✅
     - Excretion
     - `chr1:97540298-98389615 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A97540298%2D98389615&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr1:97074742-97924034 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A97074742%2D97924034&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - DRD2
     - ✅
     -
     -
     -
     -
     - Toxicity
     - `chr11:113277326-113349120 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr11%3A113277326%2D113349120&hgsid=2104698958_hLTLhdTJW6B6oZE8Chuxao0Q7W9C>`__
     - `chr11:113406604-113478398 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr11%3A113406604%2D113478398&hgsid=2104698958_hLTLhdTJW6B6oZE8Chuxao0Q7W9C>`__
     - NM_000795.4 was used as the main transcript.
   * - F2
     - ✅
     -
     -
     -
     -
     - Other
     - `chr11:46737762-46764056 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr11%3A46737762%2D46764056&hgsid=2104698958_hLTLhdTJW6B6oZE8Chuxao0Q7W9C>`__
     - `chr11:46716212-46742506 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr11%3A46716212%2D46742506&hgsid=2104698958_hLTLhdTJW6B6oZE8Chuxao0Q7W9C>`__
     - NM_000506.5 was used as the main transcript.
   * - :ref:`genes:F5`
     - ✅
     -
     - ✅
     -
     -
     - Other
     - `chr1:169478188-169558719 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A169478188%2D169558719&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr1:169508950-169589481 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A169508950%2D169589481&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:G6PD`
     - ✅
     - ✅
     - ✅
     -
     -
     - Disease
     - `chrX:153756604-153778233 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chrX%3A153756604%2D153778233&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chrX:154528389-154550018 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chrX%3A154528389%2D154550018&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - G6PD is located on X chromosome.
   * - GRIK1
     - ✅
     -
     -
     -
     -
     - Metabolism
     - `chr21:30906253-31315314 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr21%3A30906253%2D31315314&hgsid=2104698958_hLTLhdTJW6B6oZE8Chuxao0Q7W9C>`__
     - `chr21:29533932-29942996 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr21%3A29533932%2D29942996&hgsid=2104698958_hLTLhdTJW6B6oZE8Chuxao0Q7W9C>`__
     - NM_001330994.2 was used as the main transcript.
   * - GRIK4
     - ✅
     -
     -
     -
     -
     - Efficacy
     - `chr11:120379456-120862615 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr11%3A120379456%2D120862615&hgsid=2104698958_hLTLhdTJW6B6oZE8Chuxao0Q7W9C>`__
     - `chr11:120508747-120991906 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr11%3A120508747%2D120991906&hgsid=2104698958_hLTLhdTJW6B6oZE8Chuxao0Q7W9C>`__
     - NM_014619.5 was used as the main transcript.
   * - GRIN2B
     - ✅
     -
     -
     -
     -
     - Other
     - `chr12:13687270-14137536 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr12%3A13687270%2D14137536&hgsid=2115857682_DJfHD4PWNV9r1eAAvajhTGuJCkl7>`__
     - `chr12:13534336-13984602 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr12%3A13534336%2D13984602&hgsid=2115857682_DJfHD4PWNV9r1eAAvajhTGuJCkl7>`__
     - NM_000834.5 was used as the main transcript.
   * - :ref:`genes:GSTM1`
     - ✅
     - ✅
     -
     -
     -
     - Metabolism
     - `chr1:110227417-110239367 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A110227417%2D110239367&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr1:109684816-109696745 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A109684816%2D109696745&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - GSTP1
     - ✅
     -
     -
     -
     -
     - Metabolism
     - `chr11:67348065-67357124 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr11%3A67348065%2D67357124&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr11:67580811-67589653 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr11%3A67580811%2D67589653&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:GSTT1`
     -
     - ✅
     -
     -
     -
     - Metabolism
     - `chr22:24373132-24387311 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr22%3A24373132%2D24387311&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr22_KI270879v1_alt:267307-281486 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr22_KI270879v1_alt%3A267307%2D281486&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - GSTT1 is located on different contigs between GRCh37 and GRCh38.
   * - HTR1A
     - ✅
     -
     -
     -
     -
     - Efficacy
     - `chr5:63250700-63261272 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr5%3A63250700%2D63261272&hgsid=2115857682_DJfHD4PWNV9r1eAAvajhTGuJCkl7>`__
     - `chr5:63954873-63965445 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr5%3A63954873%2D63965445&hgsid=2115857682_DJfHD4PWNV9r1eAAvajhTGuJCkl7>`__
     - NM_000524.4 was used as the main transcript.
   * - HTR2A
     - ✅
     -
     -
     -
     -
     - Efficacy
     - `chr13:47402680-47474188 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr13%3A47402680%2D47474188&hgsid=2115857682_DJfHD4PWNV9r1eAAvajhTGuJCkl7>`__
     - `chr13:46828545-46900053 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr13%3A46828545%2D46900053&hgsid=2115857682_DJfHD4PWNV9r1eAAvajhTGuJCkl7>`__
     - NM_000524.4 was used as the main transcript.
   * - :ref:`genes:IFNL3`
     - ✅
     -
     - ✅
     -
     -
     - Other
     - `chr19:39731245-39744165 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A39731245%2D39744165&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr19:39240552-39253525 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A39240552%2D39253525&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - IFNL4
     - ✅
     -
     - 
     -
     -
     - Other
     - `chr19:39733953-39742496 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A39733953%2D39742496&hgsid=2144291390_4bmNOa6wq7Mk8crO6QNojSla7rfr>`__
     - `chr19:39243313-39251856 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A39243313%2D39251856&hgsid=2144291444_9qcI5Mpt1r7ap1R1ceKlp0aFJjS3>`__
     - NR_074079.1 was used as the main transcript.
   * - ITGB3
     - ✅
     -
     - 
     -
     -
     - Efficacy
     - `chr17:45328192-45394109 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr17%3A45328192%2D45394109&hgsid=2144291390_4bmNOa6wq7Mk8crO6QNojSla7rfr>`__
     - `chr17:47250826-47316743 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr17%3A47250826%2D47316743&hgsid=2144291444_9qcI5Mpt1r7ap1R1ceKlp0aFJjS3>`__
     - NM_000212.3 was used as the main transcript.
   * - ITPA
     - ✅
     -
     - 
     -
     -
     - Other
     - `chr20:3187143-3207506 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr20%3A3187143%2D3207506&hgsid=2144291390_4bmNOa6wq7Mk8crO6QNojSla7rfr>`__
     - `chr20:3206497-3226860 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr20%3A3206497%2D3226860&hgsid=2144291444_9qcI5Mpt1r7ap1R1ceKlp0aFJjS3>`__
     - NM_033453.4 was used as the main transcript.
   * - :ref:`genes:MT-RNR1`
     - ✅
     -
     - ✅
     -
     -
     - Disease
     - `chrM:410-1840 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chrMT%3A410%2D1840&hgsid=3332880867_aVXAZCwJTtLSGNhbga0xleh3qxVY>`__
     - `chrM:410-1840 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chrM%3A410%2D1840&hgsid=3332880867_aVXAZCwJTtLSGNhbga0xleh3qxVY>`__
     -
   * - MTHFR
     - ✅
     -
     -
     -
     -
     - Metabolism
     - `chr1:11842779-11869021 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A11842779%2D11869021&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr1:11782722-11808964 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A11782722%2D11808964&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - NAT1
     - ✅
     -
     -
     -
     -
     - Metabolism
     - `chr8:18064617-18084198 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr8%3A18064617%2D18084198&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr8:18207108-18226689 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr8%3A18207108%2D18226689&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - NAT2
     - ✅
     -
     -
     -
     -
     - Metabolism
     - `chr8:18245791-18261728 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr8%3A18245791%2D18261728&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr8:18388281-18404218 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr8%3A18388281%2D18404218&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:NUDT15`
     - ✅
     -
     - ✅
     - ✅
     - ✅
     - Metabolism
     - `chr13:48608702-48624364 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr13%3A48608702%2D48624364&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr13:48034725-48050221 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr13%3A48034725%2D48050221&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - OPRK1
     - ✅
     -
     -
     -
     -
     - Other
     - `chr8:54135283-54167197 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr8%3A54135283%2D54167197&hgsid=2144291390_4bmNOa6wq7Mk8crO6QNojSla7rfr>`__
     - `chr8:53222723-53254637 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr8%3A53222723%2D53254637&hgsid=2144291444_9qcI5Mpt1r7ap1R1ceKlp0aFJjS3Y>`__
     - NM_000912.5 was used as the main transcript.
   * - OPRM1
     - ✅
     -
     -
     -
     -
     - Other
     - `chr6:154357374-154456491 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr6%3A154357374%2D154456491&hgsid=2144291390_4bmNOa6wq7Mk8crO6QNojSla7rfr>`__
     - `chr6:154036239-154135356 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr6%3A154036239%2D154135356&hgsid=2144291444_9qcI5Mpt1r7ap1R1ceKlp0aFJjS3>`__
     - NM_000914.5 was used as the main transcript.
   * - :ref:`genes:POR`
     - ✅
     -
     -
     - ✅
     -
     - Disease
     - `chr7:75541419-75619173 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A75541419%2D75619173&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr7:75912154-75989855 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A75912154%2D75989855&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:PTGIS`
     - ✅
     -
     -
     - ✅
     -
     - Other
     - `chr20:48117410-48187674 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr20%3A48117410%2D48187674&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr20:49500873-49571137 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr20%3A49500873%2D49571137&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - RARG
     - ✅
     -
     - 
     - 
     -
     - Toxicity
     - `chr12:53601352-53628993 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr12%3A53601352%2D53628993&hgsid=2144291390_4bmNOa6wq7Mk8crO6QNojSla7rfr>`__
     - `chr12:53207568-53235209 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr12%3A53207568%2D53235209&hgsid=2144291444_9qcI5Mpt1r7ap1R1ceKlp0aFJjS3>`__
     - NM_000966.6 was used as the main transcript.
   * - :ref:`genes:RYR1`
     - ✅
     -
     - ✅
     - ✅
     -
     - Disease
     - `chr19:38921339-39081204 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A38921339%2D39081204&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr19:38430690-38590564 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A38430690%2D38590564&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - SLC6A4
     - ✅
     -
     -
     -
     -
     - Other
     - `chr17:28518336-28565715 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr17%3A28518336%2D28565715&hgsid=2144291390_4bmNOa6wq7Mk8crO6QNojSla7rfr>`__
     - `chr17:30191318-30238697 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr17%3A30191318%2D30238697&hgsid=2144291444_9qcI5Mpt1r7ap1R1ceKlp0aFJjS3>`__
     - NM_001045.6 was used as the main transcript.
   * - SLC15A2
     - ✅
     -
     -
     -
     -
     - Excretion
     - `chr3:121610170-121666034 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr3%3A121610170%2D121666034&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr3:121891400-121947188 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr3%3A121891400%2D121947188&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:SLC22A2`
     - ✅
     - ✅
     -
     -
     -
     - Excretion
     - `chr6:160627786-160689853 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr6%3A160627786%2D160689853&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr6:160206754-160268821 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr6%3A160206754%2D160268821&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - SLC28A3
     - ✅
     - ✅
     -
     -
     -
     - Excretion
     - `chr9:86887371-86958673 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr9%3A86887371%2D86958673&hgsid=2144291390_4bmNOa6wq7Mk8crO6QNojSla7rfr>`__
     - `chr9:84272456-84343758 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr9%3A84272456%2D84343758&hgsid=2144291444_9qcI5Mpt1r7ap1R1ceKlp0aFJjS3>`__
     - NM_001199633.2 was used as the main transcript.
   * - SLC47A2
     - ✅
     - 
     -
     -
     -
     - Other
     - `chr17:19578629-19622897 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr17%3A19578629%2D19622897&hgsid=2144291390_4bmNOa6wq7Mk8crO6QNojSla7rfr>`__
     - `chr17:19675316-19719584 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr17%3A19675316%2D19719584&hgsid=2144291444_9qcI5Mpt1r7ap1R1ceKlp0aFJjS3>`__
     - NM_001099646.3 was used as the main transcript.
   * - :ref:`genes:SLCO1B1`
     - ✅
     -
     - ✅
     - ✅
     - ✅
     - Absorption
     - `chr12:21281127-21395730 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr12%3A21281127%2D21395730&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr12:21128193-21242796 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr12%3A21128193%2D21242796&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - SLCO1B3
     - ✅
     -
     -
     -
     -
     - Absorption
     - `chr12:20960637-21072845 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr12%3A20960637%2D21072845&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr12:20807704-20919911 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr12%3A20807704%2D20919911&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - SLCO2B1
     - ✅
     -
     -
     -
     -
     - Absorption
     - `chr11:74859151-74920594 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr11%3A74859151%2D74920594&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr11:75148106-75209549 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr11%3A75148106%2D75209549&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:SULT1A1`
     - ✅
     - ✅
     -
     -
     -
     - Metabolism
     - `chr16:28601907-28636365 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr16%3A28601907%2D28636365&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr16:28590586-28625044 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr16%3A28590586%2D28625044&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:TBXAS1`
     - ✅
     -
     -
     - ✅
     -
     - Other
     - `chr7:139525951-139723125 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A139525951%2D139723125&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr7:139826263-140023321 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A139826263%2D140023321&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:TPMT`
     - ✅
     -
     - ✅
     -
     - ✅
     - Metabolism
     - `chr6:18125541-18158400 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr6%3A18125541%2D18158400&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr6:18125310-18158169 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr6%3A18125310%2D18158169&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:UGT1A1`
     - ✅
     -
     - ✅
     -
     - ✅
     - Excretion
     - `chr2:234662918-234687945 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2%3A234662918%2D234687945&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr2:233754269-233779300 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2%3A233754269%2D233779300&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:UGT1A4`
     - ✅
     - ✅
     -
     -
     -
     - Excretion
     - `chr2:234624437-234684945 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2%3A234624437%2D234684945&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr2:233715735-233776300 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2%3A233715735%2D233776300&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - UGT1A6
     - ✅
     - 
     -
     -
     -
     - Excretion
     - `chr2:234598535-234684946 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2%3A234598535%2D234684946&hgsid=2144291390_4bmNOa6wq7Mk8crO6QNojSla7rfr>`__
     - `chr2:233689889-233776300 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2%3A233689889%2D233776300&hgsid=2144291444_9qcI5Mpt1r7ap1R1ceKlp0aFJjS3>`__
     - NM_001072.4 was used as the main transcript.
   * - UGT2B7
     - ✅
     -
     -
     -
     -
     - Excretion
     - `chr4:69959191-69981705 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr4%3A69959191%2D69981705&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr4:69093473-69115987 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr4%3A69093473%2D69115987&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:UGT2B15`
     - ✅
     - ✅
     -
     -
     -
     - Excretion
     - `chr4:69506314-69542494 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr4%3A69506314%2D69542494&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr4:68640596-68676652 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr4%3A68640596%2D68676652&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:UGT2B17`
     -
     - ✅
     -
     -
     -
     - Excretion
     - `chr4:69399901-69437245 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr4%3A69399901%2D69437245&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr4:68534183-68571527 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr4%3A68534183%2D68571527&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - VKORC1
     - ✅
     -
     -
     -
     - ✅
     - Target
     - `chr16:31099162-31109320 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr16%3A31099162%2D31109320&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr16:31087853-31097797 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr16%3A31087853%2D31097797&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - XPC
     - ✅
     -
     -
     -
     -
     - Other
     - `chr3:14183646-14223172 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr3%3A14183646%2D14223172&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr3:14142146-14181672 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr3%3A14142146%2D14181672&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -

ABCG2
=====

Phenotype summary for ABCG2
---------------------------

Diplotype-phenotype mapping is used for phenotype prediction.

.. list-table::
   :header-rows: 1

   * - Phenotype
     - Example
     - Priority
   * - Normal Function
     - Reference/Reference
     - Normal/Routine/Low Risk
   * - Decreased Function
     - Reference/rs2231142
     - Abnormal/Priority/High Risk
   * - Poor Function
     - rs2231142/rs2231142
     - Abnormal/Priority/High Risk

Recommendations for ABCG2
-------------------------

- Rosuvastatin

  "Prescribe ≤20mg as a starting dose and adjust doses of rosuvastatin based
  on disease-specific and specific population guidelines for patients who are
  SLCO1B1 or ABCG2 poor function phenotype. If dose >20mg needed for desired
  efficacy, consider combination therapy (i.e. rosuvastatin plus non-statin
  guideline directed medical therapy). Patients with both ABCG2 poor function
  and SLCO1B1 poor/decreased function should be prescribed ≤10mg as a
  starting dose." (Source: `PharmGKB <https://www.pharmgkb.org/
  guidelineAnnotation/PA166262321>`__)

Resources for ABCG2
-------------------

- `CPIC® guideline for statins and SLCO1B1, ABCG2, and CYP2C9 <https://cpicpgx.org/guidelines/cpic-guideline-for-statins/>`__
- `The Clinical Pharmacogenetics Implementation Consortium Guideline for SLCO1B1, ABCG2, and CYP2C9 genotypes and Statin-Associated Musculoskeletal Symptoms <https://ascpt.onlinelibrary.wiley.com/doi/10.1002/cpt.2557>`__
- `PharmGKB: Annotation of CPIC Guideline for rosuvastatin and ABCG2, SLCO1B1 <https://www.pharmgkb.org/guidelineAnnotation/PA166262321>`__

CACNA1S
=======

Phenotype summary for CACNA1S
-----------------------------

Diplotype-phenotype mapping is used for phenotype prediction.

 .. list-table::
    :header-rows: 1

    * - Phenotype
      - Example
      - Priority
    * - Uncertain Susceptibility
      - Reference/Reference
      - Normal Risk
    * - Malignant Hyperthermia Susceptibility
      - Reference/c.520C>T
      - Abnormal/Priority/High Risk

Recommendations for CACNA1S
---------------------------

- Desflurane

  "The CPIC Dosing Guideline recommends that halogenated volatile anesthetics
  uch as desflurane, enflurane, halothane, isoflurane, methoxyflurane,
  sevoflurane and the depolarizing muscle relaxants succinylcholine are
  relatively contraindicated in persons with malignant hyperthermia
  susceptibility (MHS). See full guideline for disclaimers, further details
  and supporting evidence." (Source: `PharmGKB <https://www.pharmgkb.org/
  chemical/PA164749136/guidelineAnnotation/PA166180457>`__)

- Enflurane

  "The CPIC Dosing Guideline recommends that halogenated volatile anesthetics
  such as desflurane, enflurane, halothane, isoflurane, methoxyflurane,
  sevoflurane and the depolarizing muscle relaxants succinylcholine are
  relatively contraindicated in persons with malignant hyperthermia
  susceptibility (MHS). See full guideline for disclaimers, further details
  and supporting evidence." (Source: `PharmGKB <https://www.pharmgkb.org/
  chemical/PA449461/guidelineAnnotation/PA166180457>`__)

- Halothane

  "The CPIC Dosing Guideline recommends that halogenated volatile anesthetics
  such as desflurane, enflurane, halothane, isoflurane, methoxyflurane,
  sevoflurane and the depolarizing muscle relaxants succinylcholine are
  relatively contraindicated in persons with malignant hyperthermia
  susceptibility (MHS). See full guideline for disclaimers, further details
  and supporting evidence." (Source: `PharmGKB <https://www.pharmgkb.org/
  chemical/PA449845/guidelineAnnotation/PA166180457>`__)

- Isoflurane

  "The CPIC Dosing Guideline recommends that halogenated volatile anesthetics
  such as desflurane, enflurane, halothane, isoflurane, methoxyflurane,
  sevoflurane and the depolarizing muscle relaxants succinylcholine are
  relatively contraindicated in persons with malignant hyperthermia
  susceptibility (MHS). See full guideline for disclaimers, further details
  and supporting evidence." (Source: `PharmGKB <https://www.pharmgkb.org/
  chemical/PA450106/guidelineAnnotation/PA166180457>`__)

- Methoxyflurane

  "The CPIC Dosing Guideline recommends that halogenated volatile anesthetics
  such as desflurane, enflurane, halothane, isoflurane, methoxyflurane,
  sevoflurane and the depolarizing muscle relaxants succinylcholine are
  relatively contraindicated in persons with malignant hyperthermia
  susceptibility (MHS). See full guideline for disclaimers, further details
  and supporting evidence." (Source: `PharmGKB <https://www.pharmgkb.org/
  chemical/PA450434/guidelineAnnotation/PA166180457>`__)

- Sevoflurane

  "The CPIC Dosing Guideline recommends that halogenated volatile anesthetics
  such as desflurane, enflurane, halothane, isoflurane, methoxyflurane,
  sevoflurane and the depolarizing muscle relaxants succinylcholine are
  relatively contraindicated in persons with malignant hyperthermia
  susceptibility (MHS). See full guideline for disclaimers, further details
  and supporting evidence." (Source: `PharmGKB <https://www.pharmgkb.org/
  chemical/PA451341/guidelineAnnotation/PA166180457>`__)

- Succinylcholine

  "The CPIC Dosing Guideline recommends that halogenated volatile anesthetics
  such as desflurane, enflurane, halothane, isoflurane, methoxyflurane,
  sevoflurane and the depolarizing muscle relaxants succinylcholine are
  relatively contraindicated in persons with malignant hyperthermia
  susceptibility (MHS). See full guideline for disclaimers, further details
  and supporting evidence." (Source: `PharmGKB <https://www.pharmgkb.org/
  chemical/PA451522/guidelineAnnotation/PA166180457>`__)

Resources for CACNA1S
---------------------

- `PharmGKB: Annotation of CPIC Guideline for desflurane and CACNA1S, RYR1 <https://www.pharmgkb.org/chemical/PA164749136/guidelineAnnotation/PA166180457>`__
- `CPIC® Guideline for Potent Volatile Anesthetic Agents and Succinylcholine and RYR1 and CACNA1S <https://cpicpgx.org/guidelines/cpic-guideline-for-ryr1-and-cacna1s/>`__
- `Clinical Pharmacogenetics Implementation Consortium (CPIC) Guideline for the Use of Potent Volatile Anesthetic Agents and Succinylcholine in the Context of RYR1 or CACNA1S Genotypes <https://doi.org/10.1002/cpt.1319>`__

CFTR
====

Phenotype summary for CFTR
--------------------------

Diplotype-phenotype mapping is used for phenotype prediction.

 .. list-table::
    :header-rows: 1

    * - Phenotype
      - Example
      - Priority
    * - Favorable Response
      - Reference/G551D
      - None
    * - Unfavorable Response
      - F508del/F508del
      - None
    * - Indeterminate
      - Reference/F508del
      - None

Recommendations for CFTR
------------------------

- Ivacaftor

  "Ivacaftor treatment is recommended only in cystic fibrosis (CF) patients
  that are either homozygous or heterozygous for certain CFTR variants. See
  full guideline for disclaimers, further details and supporting evidence."
  (Source: `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/
  PA166114461>`__)

Resources for CFTR
------------------

- `PharmGKB: Annotation of CPIC Guideline for ivacaftor and CFTR <https://www.pharmgkb.org/chemical/PA165950341/guidelineAnnotation/PA166114461>`__
- `CPIC® Guideline for Ivacaftor and CFTR <https://cpicpgx.org/guidelines/guideline-for-ivacaftor-and-cftr/>`__
- `Clinical Pharmacogenetics Implementation Consortium (CPIC) Guidelines for Ivacaftor Therapy in the Context of CFTR Genotype <https://doi.org/10.1038/clpt.2014.54>`__

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

   * - SV Alleles
     - SV Name
     - Genotype
     - Reference
     - Gene Model
     - GRCh37
     - GRCh38
     - Data Type
     - Source
     - Coriell ID
     - Version
     - Description
   * -
     - Normal
     - \*1/\*2
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-1.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-5.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-5.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA10831
     - 0.4.1
     -
   * - \*4
     - WholeDel1
     - \*1/\*4
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-2.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-1.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-1.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA18617
     - 0.4.1
     -
   * - \*4
     - WholeDel1Hom
     - \*4/\*4
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-3.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-2.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-2.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA18952
     - 0.4.1
     -
   * - \*4
     - WholeDel2
     - \*1/\*4
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-2.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-6.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-6.png>`
     - WGS
     -
     -
     - 0.12.0
     -
   * - \*4
     - WholeDel2Hom
     - \*4/\*4
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-3.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-16.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-16.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA21093
     - 0.15.0
     -
   * - \*4
     - WholeDel3
     - \*4/\*9
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-2.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-7.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-7.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA18488
     - 0.12.0
     -
   * - \*1x2
     - WholeDup1
     - \*1x2/\*25
     - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-4.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-3.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-3.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA18861
     - 0.4.1
     -
   * - \*1x2
     - WholeDup2
     - \*1x2/\*2
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-4.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-10.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-10.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA12342
     - 0.12.0
     -
   * - \*1x2
     - WholeDup3
     - \*1x2/\*17
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-4.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-11.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-11.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA19129
     - 0.12.0
     -
   * -
     - Hybrid1
     - Indeterminate
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-11.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-4.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-4.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - HG00436
     - 0.4.1
     -
   * - \*12
     - Hybrid2
     - \*1/\*12
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-5.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-8.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-8.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA11829
     - 0.12.0
     - \*12 has exons 1-2 of CYP2A7 origin and exons 3-9 of CYP2A6 origin (breakpoint in intron 2).
   * - \*12
     - Hybrid2Hom
     - \*12/\*12
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-9.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-14.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-14.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA19780
     - 0.14.0
     -
   * - \*34
     - Hybrid3
     - \*1/\*34
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-6.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-9.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-9.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA18516
     - 0.12.0
     - \*34 has exons 1-4 of CYP2A7 origin and exons 5-9 of CYP2A6 origin (breakpoint in intron 4).
   * -
     - Hybrid4
     - Indeterminate
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-10.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-15.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-15.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA20515
     - 0.14.0
     -
   * -
     - Hybrid5
     - Indeterminate
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-13.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-17.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-17.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - HG00155
     - 0.15.0
     -
   * -
     - Hybrid6
     - Indeterminate
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-12.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-18.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-18.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - HG00141
     - 0.15.0
     -
   * -
     - Hybrid7
     - Indeterminate
     -
     -
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-21.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-21.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - HG02382
     - 0.16.0
     -
   * -
     - Tandem1
     - Indeterminate
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-8.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-13.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-13.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA20828
     - 0.14.0
     -
   * -
     - Tandem2
     - Indeterminate
     -
     -
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-22.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-22.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - HG04214
     - 0.16.0
     -
   * -
     - ParalogWholeDel1
     - Indeterminate
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-14.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-19.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-19.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - HG00625
     - 0.15.0
     -
   * -
     - ParalogWholeDup1
     - Indeterminate
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-7.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-12.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-12.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA06985
     - 0.12.0
     -
   * -
     - Unknown1
     - Indeterminate
     -
     -
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-20.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-20.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - HG02081
     - 0.16.0
     -

PyPGx was recently applied to the entire high-coverage WGS dataset from 1KGP
(N=2,504). Click `here <https://github.com/sbslee/1kgp-pgx-paper/blob/main/
sv-tables/CYP2A6.md>`__ to see individual SV calls for CYP2A6, and
corresponding copy number profiles and allele fraction profiles.

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

Resources for CYP2A6
--------------------

- `PharmVar CYP2A6 page <https://www.pharmvar.org/gene/CYP2A6>`__

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

   * - SV Alleles
     - SV Name
     - Genotype
     - Reference
     - Gene Model
     - GRCh37
     - GRCh38
     - Data Type
     - Source
     - Coriell ID
     - Version
     - Description
   * -
     - Normal
     - \*1/\*2
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2B6-1.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2B6-2.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2B6-2.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA12813
     - 0.4.1
     -
   * - \*22x2
     - WholeDup1
     - \*6/\*22x2
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2B6-3.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2B6-3.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2B6-3.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA19190
     - 0.12.0
     -
   * - \*29
     - Hybrid1
     - \*6/\*29
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2B6-2.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2B6-1.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2B6-1.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA19178
     - 0.4.1
     - \*29 has exons 1-4 of CYP2B7 origin and exons 5-9 of CYP2A6 origin (breakpoint in intron 4).
   * -
     - Tandem1
     - Indeterminate
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2B6-4.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2B6-4.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2B6-4.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - HG01806
     - 0.16.0
     -
   * -
     - PartialDup1
     - Indeterminate
     -
     -
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2B6-5.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2B6-5.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - HG03784
     - 0.16.0
     -
   * -
     - PartialDup2
     - Indeterminate
     -
     -
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2B6-6.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2B6-6.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - HG02790
     - 0.16.0
     -
   * -
     - ParalogWholeDel1
     - Indeterminate
     -
     -
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2B6-7.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2B6-7.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - HG03235
     - 0.16.0
     -

PyPGx was recently applied to the entire high-coverage WGS dataset from 1KGP
(N=2,504). Click `here <https://github.com/sbslee/1kgp-pgx-paper/blob/main/
sv-tables/CYP2B6.md>`__ to see individual SV calls for CYP2B6, and
corresponding copy number profiles and allele fraction profiles.

Phenotype summary for CYP2B6
----------------------------

Diplotype-phenotype mapping is used for phenotype prediction.

 .. list-table::
    :header-rows: 1

    * - Phenotype
      - Example
      - Priority
    * - Ultrarapid Metabolizer
      - \*4/\*4
      - Normal/Routine/Low Risk
    * - Rapid Metabolizer
      - \*1/\*4
      - Normal/Routine/Low Risk
    * - Normal Metabolizer
      - \*1/\*2
      - Normal/Routine/Low Risk
    * - Intermediate Metabolizer
      - \*1/\*29
      - Abnormal/Priority/High Risk
    * - Poor Metabolizer
      - \*6/\*6
      - Abnormal/Priority/High Risk
    * - Indeterminate
      - \*1/\*3
      - None

Recommendations for CYP2B6
--------------------------

- Efavirenz

  "Consider initiating efavirenz with a decreased dose of either 400 or 200
  mg/day for patients who are CYP2B6 poor metabolizers. Consider initiating
  efavirenz with a decreased dose of 400 mg/day for patients who are CYP2B6
  intermediate metabolizers." (Source: `PharmGKB <https://www.pharmgkb.org
  /guidelineAnnotation/PA166182603>`__)

Resources for CYP2B6
--------------------

- `PharmVar CYP2B6 page <https://www.pharmvar.org/gene/CYP2B6>`__
- `CPIC® Guideline for Efavirenz based on CYP2B6 genotype <https://cpicpgx.org/guidelines/cpic-guideline-for-efavirenz-based-on-cyp2b6-genotype/>`__
- `PharmGKB: Annotation of CPIC Guideline for efavirenz and CYP2B6 <https://www.pharmgkb.org/guidelineAnnotation/PA166182603>`__

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
     - Priority
   * - Normal Metabolizer
     - 2 == score
     - \*1/\*1
     - Normal/Routine/Low Risk
   * - Intermediate Metabolizer
     - 1 <= score < 2
     - \*1/\*2
     - Abnormal/Priority/High Risk
   * - Poor Metabolizer
     - 0 <= score < 1
     - \*2/\*3
     - Abnormal/Priority/High Risk
   * - Indeterminate
     - Unknown score
     - \*1/\*7
     - Normal/Routine/Low Risk

Recommendations for CYP2C9
--------------------------

- Celecoxib

  "The CPIC Dosing Guideline for celecoxib, flurbiprofen, ibuprofen and
  lornoxicam recommends initiating therapy with 25-50% of the lowest
  recommended starting dose for CYP2C9 poor metabolizers and initiating
  therapy with lowest recommended starting dose for CYP2C9 intermediate
  metabolizers with activity score of 1. See full guideline for further
  details and supporting evidence." (Source: `PharmGKB <https://www.pharmgkb.
  org/chemical/PA448871/guidelineAnnotation/PA166191841>`__)

- Flurbiprofen

  "The CPIC Dosing Guideline for celecoxib, flurbiprofen, ibuprofen and
  lornoxicam recommends initiating therapy with 25-50% of the lowest
  recommended starting dose for CYP2C9 poor metabolizers and initiating
  therapy with lowest recommended starting dose for CYP2C9 intermediate
  metabolizers with activity score of 1. See full guideline for further
  details and supporting evidence." (Source: `PharmGKB <https://www.pharmgkb.
  org/chemical/PA449683/guidelineAnnotation/PA166191841>`__)

- Fluvastatin

  "CYP2C9 IMs should avoid fluvastatin doses greater than 40mg while CYP2C9
  PMs should avoid doses greater than 20mg. If higher doses are required for
  desired efficacy, an alternative statin should be considered. Patients with
  SLCO1B1 poor function should also avoid fluvastatin doses greater than 40mg
  and and adjust doses of fluvastatin based on disease-specific guidelines.
  Patients with both SLCO1B1 poor function and CYP2C9 IM/PM should be
  prescribed an alternative statin depending on the desired potency."
  (Source: `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/
  PA166262341>`__)

- Ibuprofen

  "The CPIC Dosing Guideline for celecoxib, flurbiprofen, ibuprofen and
  lornoxicam recommends initiating therapy with 25-50% of the lowest
  recommended starting dose for CYP2C9 poor metabolizers and initiating
  therapy with lowest recommended starting dose for CYP2C9 intermediate
  metabolizers with activity score of 1. See full guideline for further
  details and supporting evidence." (Source: `PharmGKB <https://www.pharmgkb.
  org/chemical/PA449957/guidelineAnnotation/PA166191841>`__)

- Lornoxicam

  "The CPIC Dosing Guideline for celecoxib, flurbiprofen, ibuprofen and
  lornoxicam recommends initiating therapy with 25-50% of the lowest
  recommended starting dose for CYP2C9 poor metabolizers and initiating
  therapy with lowest recommended starting dose for CYP2C9 intermediate
  metabolizers with activity score of 1. See full guideline for further
  details and supporting evidence." (Source: `PharmGKB <https://www.pharmgkb.
  org/chemical/PA165958395/guidelineAnnotation/PA166191841>`__)

- Meloxicam

  "The CPIC Dosing Guideline for meloxicam recommends alternative therapy for
  CYP2C9 poor metabolizers due to markedly prolonged half-life, and
  initiating therapy with 50% of the lowest recommended starting dose or
  choose an alternative therapy for CYP2C9 intermediate metabolizers with
  activity score of 1. See full guideline for further details and supporting
  evidence." (Source: `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation
  /PA166192301>`__)

- Piroxicam

  "The CPIC Dosing Guideline for piroxicam recommends that CYP2C9 poor
  metabolizers and intermediate metabolizers with activity score of 1 should
  choose an alternative therapy not metabolized by CYP2C9 or not
  significantly impacted by CYP2C9 genetic variants in vivo or choose an
  NSAID metabolized by CYP2C9 but with a shorter half-life. See full
  guideline for further details and supporting evidence." (Source:
  `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/PA166192321>`__)

- Tenoxicam

  "The CPIC Dosing Guideline for tenoxicam recommends that CYP2C9 poor
  metabolizers and intermediate metabolizers with activity score of 1 should
  choose an alternative therapy not metabolized by CYP2C9 or not
  significantly impacted by CYP2C9 genetic variants in vivo or choose an
  NSAID metabolized by CYP2C9 but with a shorter half-life. See full
  guideline for further details and supporting evidence." (Source:
  `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/PA166192341>`__)

Resources for CYP2C9
--------------------

- `PharmVar CYP2C9 page <https://www.pharmvar.org/gene/CYP2C9>`__
- `CPIC® Guideline for NSAIDs based on CYP2C9 genotype <https://cpicpgx.org/guidelines/cpic-guideline-for-nsaids-based-on-cyp2c9-genotype/>`__
- `The Clinical Pharmacogenetics Implementation Consortium Guideline for SLCO1B1, ABCG2, and CYP2C9 genotypes and Statin-Associated Musculoskeletal Symptoms <https://ascpt.onlinelibrary.wiley.com/doi/10.1002/cpt.2557>`__

CYP2C19
=======

Phenotype summary for CYP2C19
-----------------------------

Diplotype-phenotype mapping is used for phenotype prediction.

.. list-table::
   :header-rows: 1

   * - Phenotype
     - Example
     - Priority
   * - Ultrarapid Metabolizer
     - \*17/\*17
     - Abnormal/Priority/High Risk
   * - Rapid Metabolizer
     - \*1/\*17
     - Abnormal/Priority/High Risk
   * - Normal Metabolizer
     - \*1/\*1
     - Normal/Routine/Low Risk
   * - Likely Intermediate Metabolizer
     - \*1/\*10
     - Abnormal/Priority/High Risk
   * - Intermediate Metabolizer
     - \*1/\*2
     - Abnormal/Priority/High Risk
   * - Likely Poor Metabolizer
     - \*10/\*22
     - Abnormal/Priority/High Risk
   * - Poor Metabolizer
     - \*2/\*2
     - Abnormal/Priority/High Risk
   * - Indeterminate
     - \*1/\*12
     - None

Recommendations for CYP2C19
---------------------------

- Amitriptyline

  "The CPIC Dosing Guideline update for amitriptyline recommends an
  alternative drug for CYP2D6 ultrarapid or poor metabolizers and CYP2C19
  ultrarapid, rapid or poor metabolizers. If amitriptyline is warranted,
  consider a 50% dose reduction in CYP2D6 or CYP2C19 poor metabolizers. For
  CYP2D6 intermediate metabolizers, a 25% dose reduction should be
  considered." (Source: `PharmGKB <https://www.pharmgkb.org
  /guidelineAnnotation/PA166105006>`__)

- Citalopram

  "The CPIC Dosing Guideline for the selective serotonin reuptake inhibitors
  citalopram and escitalopram recommends an alternative drug not
  predominantly metabolized by CYP2C19 for CYP2C19 ultrarapid metabolizers.
  For CYP2C19 poor metabolizers, consider a 50% reduction of recommended
  starting dose and titrate to response or select alternative drug not
  predominantly metabolized by CYP2C19." (Source: `PharmGKB
  <https://www.pharmgkb.org/chemical/PA449015/guidelineAnnotation/
  PA166127638>`__)

- Clomipramine

  "Tricyclic antidepressants have comparable pharmacokinetic properties, it
  may be reasonable to apply the CPIC Dosing Guideline for amitriptyline and
  CYP2C19, CYP2D6 to other tricyclics including clomipramine. The CPIC Dosing
  Guideline update for amitriptyline recommends an alternative drug for
  CYP2D6 ultrarapid or poor metabolizers and CYP2C19 ultrarapid, rapid or
  poor metabolizers. If amitriptyline is warranted, consider a 50% dose
  reduction in CYP2D6 or CYP2C19 poor metabolizers. For CYP2D6 intermediate
  metabolizers, a 25% dose reduction should be considered." (Source:
  `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/PA166105007>`__)

- Clopidogrel

  "The CPIC Dosing Guideline for clopidogrel recommends an alternative
  antiplatelet therapy for CYP2C19 poor or intermediate metabolizers
  (cardiovascular indications: prasugrel or ticagrelor if no
  contraindication; neurovascular indications: alternative P2Y12 inhibitor
  if clinically indicated and no contraindication.)" (Source: `PharmGKB
  <https://www.pharmgkb.org/guidelineAnnotation/PA166104948>`__)

- Dexlansoprazole

  "The CPIC Dosing Guideline recommendations for dexlansoprazole are based on
  the similarity in its metabolism and lansoprazole and extrapolated from the
  first-generation PPIs. The guideline recommends to increase the starting
  daily dose and to monitor efficacy in CYP2C19 ultrarapid metabolizer. For
  CYP2C19 rapid and normal metabolizers in the treatment of H. pylori
  infection and erosive esophagitis increasing the dose might be considered
  after initiation with the standard starting daily dose. The recommendations
  for intermediate and poor metabolizers for chronic therapy (>12 weeks) and
  efficacy achieved is to consider 50% reduction in daily dose. See full
  guideline for further details and supporting evidence." (Source: `PharmGKB
  <https://www.pharmgkb.org/guidelineAnnotation/PA166219301>`__)

- Doxepin

  "Tricyclic antidepressants have comparable pharmacokinetic properties, it
  may be reasonable to apply the CPIC Dosing Guideline for amitriptyline and
  CYP2C19, CYP2D6 to other tricyclics including doxepin. The CPIC Dosing
  Guideline update for amitriptyline recommends an alternative drug for
  CYP2D6 ultrarapid or poor metabolizers and CYP2C19 ultrarapid, rapid or
  poor metabolizers. If amitriptyline is warranted, consider a 50% dose
  reduction in CYP2D6 or CYP2C19 poor metabolizers. For CYP2D6 intermediate
  metabolizers, a 25% dose reduction should be considered." (Source:
  `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/PA166105000>`__)

- Escitalopram

  "The CPIC Dosing Guideline for the selective serotonin reuptake inhibitors
  citalopram and escitalopram recommends an alternative drug not
  predominantly metabolized by CYP2C19 for CYP2C19 ultrarapid metabolizers.
  For CYP2C19 poor metabolizers, consider a 50% reduction of recommended
  starting dose and titrate to response or select alternative drug not
  predominantly metabolized by CYP2C19." (Source: `PharmGKB
  <https://www.pharmgkb.org/chemical/PA10074/guidelineAnnotation/
  PA166127638>`__)

- Imipramine

  "Tricyclic antidepressants have comparable pharmacokinetic properties, it
  may be reasonable to apply the CPIC Dosing Guideline for amitriptyline and
  CYP2C19, CYP2D6 to other tricyclics including imipramine. The CPIC Dosing
  Guideline update for amitriptyline recommends an alternative drug for
  CYP2D6 ultrarapid or poor metabolizers and CYP2C19 ultrarapid, rapid or
  poor metabolizers. If amitriptyline is warranted, consider a 50% dose
  reduction in CYP2D6 or CYP2C19 poor metabolizers. For CYP2D6 intermediate
  metabolizers, a 25% dose reduction should be considered." (Source:
  `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/PA166104999>`__)

- Lansoprazole

  "The CPIC Dosing Guideline for omeprazole, lansoprazole, pantoprazole
  recommends to increase the starting daily dose and to monitor efficacy in
  CYP2C19 ultrarapid metabolizer. For CYP2C19 rapid and normal metabolizers in
  the treatment of H. pylori infection and erosive esophagitis increasing the
  dose might be considered after initiation with the standard starting daily
  dose. The recommendations for intermediate and poor metabolizer for chronic
  therapy (>12 weeks) and efficacy achieved is to consider 50% reduction in
  daily dose. See full guideline for further details and supporting evidence."
  (Source: `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/
  PA166219103>`__)

- Omeprazole

  "The CPIC Dosing Guideline for omeprazole, lansoprazole, pantoprazole
  recommends to increase the starting daily dose and to monitor efficacy in
  CYP2C19 ultrarapid metabolizer. For CYP2C19 rapid and normal metabolizers in
  the treatment of H. pylori infection and erosive esophagitis increasing the
  dose might be considered after initiation with the standard starting daily
  dose. The recommendations for intermediate and poor metabolizer for chronic
  therapy (>12 weeks) and efficacy achieved is to consider 50% reduction in
  daily dose. See full guideline for further details and supporting evidence."
  (Source: `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/
  PA166219103>`__)

- Pantoprazole

  "The CPIC Dosing Guideline for omeprazole, lansoprazole, pantoprazole
  recommends to increase the starting daily dose and to monitor efficacy in
  CYP2C19 ultrarapid metabolizer. For CYP2C19 rapid and normal metabolizers in
  the treatment of H. pylori infection and erosive esophagitis increasing the
  dose might be considered after initiation with the standard starting daily
  dose. The recommendations for intermediate and poor metabolizer for chronic
  therapy (>12 weeks) and efficacy achieved is to consider 50% reduction in
  daily dose. See full guideline for further details and supporting evidence."
  (Source: `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/
  PA166219103>`__)

- Sertraline

  "The CPIC Dosing Guideline for the selective serotonin reuptake inhibitor
  sertraline recommends to consider a 50% reduction of recommended starting
  dose and titrate to response or select alternative drug not predominantly
  metabolized by CYP2C19 for CYP2C19 poor metabolizers." (Source: `PharmGKB
  <https://www.pharmgkb.org/guidelineAnnotation/PA166127639>`__)

- Voriconazole

  "The CPIC dosing guideline for voriconazole recommends selecting an
  alternative agent that is not dependent on CYP2C19 metabolism in adults who
  are CYP2C19 ultrarapid metabolizers, rapid metabolizers or poor
  metabolizers. In pediatric patients, an alternative agent should be used in
  patients who are ultrarapid metabolizers or poor metabolizers. In pediatric
  rapid metabolizers, therapy should be initiated at recommended standard
  case dosing, then therapeutic dosing monitoring should be used to titrate
  dose to therapeutic trough concentrations." (Source: `PharmGKB
  <https://www.pharmgkb.org/guidelineAnnotation/PA166161537>`__)

- Trimipramine

  "Tricyclic antidepressants have comparable pharmacokinetic properties, it
  may be reasonable to apply the CPIC Dosing Guideline for amitriptyline and
  CYP2C19, CYP2D6 to other tricyclics including trimipramine. The CPIC Dosing
  Guideline update for amitriptyline recommends an alternative drug for
  CYP2D6 ultrarapid or poor metabolizers and CYP2C19 ultrarapid, rapid or
  poor metabolizers. If amitriptyline is warranted, consider a 50% dose
  reduction in CYP2D6 or CYP2C19 poor metabolizers. For CYP2D6 intermediate
  metabolizers, a 25% dose reduction should be considered." (Source:
  `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/PA166105001>`__)

Resources for CYP2C19
---------------------

- `PharmVar CYP2C19 page <https://www.pharmvar.org/gene/CYP2C19>`__
- `CPIC® Guideline for Voriconazole and CYP2C19 <https://cpicpgx.org/guidelines/guideline-for-voriconazole-and-cyp2c19/>`__

CYP2D6
======

SV summary for CYP2D6
---------------------

Below is comprehensive summary of SV described from real NGS studies:

.. list-table::
   :header-rows: 1

   * - SV Alleles
     - SV Name
     - Genotype
     - Reference
     - Gene Model
     - GRCh37
     - GRCh38
     - Data Type
     - Source
     - Coriell ID
     - Version
     - Description
   * -
     - Normal
     - \*1/\*2
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2D6-1.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-8.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-8.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA11839
     - 0.4.1
     -
   * - \*5
     - WholeDel1
     - \*5/\*29
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2D6-2.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-1.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-1.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA18861
     - 0.4.1
     -
   * - \*5
     - WholeDel1Hom
     - \*5/\*5
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2D6-3.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-6.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-6.png>`
     - WGS
     -
     -
     - 0.10.0
     -
   * - \*4x2
     - WholeDup1
     - \*2/\*4x2
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2D6-4.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-2.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-2.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA19819
     - 0.4.1
     -
   * - \*1x3
     - WholeMultip1
     - \*1x3/\*10
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2D6-5.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-12.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-12.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA19190
     - 0.12.0
     -
   * - \*68+\*4
     - Tandem1A
     - \*139/\*68+\*4
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2D6-6.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-3.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-3.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA11832
     - 0.4.1
     - \*68 has exon 1 of CYP2D6 origin and exons 2-9 of CYP2D7 origin (breakpoint in intron 1).
   * - \*68+\*4
     - Tandem1B
     - \*68+\*4/\*68+\*4
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2D6-7.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-13.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-13.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA12282
     - 0.12.0
     -
   * - \*36+\*10
     - Tandem2A
     - \*2/\*36+\*10
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2D6-8.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-4.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-4.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA18564
     - 0.4.1
     - \*36 has exons 1-8 of CYP2D6 origin and exon 9 of CYP2D7 origin (breakpoint in exon 9).
   * - \*36x2+\*10
     - Tandem2B
     - \*1/\*36x2+\*10
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2D6-9.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-5.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-5.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA18524
     - 0.4.1
     -
   * - \*36x3+\*10
     - Tandem2C
     - \*1/\*36x3+\*10
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2D6-10.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-7.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-7.png>`
     - WGS
     -
     -
     - 0.10.0
     -
   * -
     - Tandem2F
     - Indeterminate
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2D6-16.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-19.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-19.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - HG00458
     - 0.15.0
     -
   * - \*13+\*1
     - Tandem3
     - \*1/\*13+\*1
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2D6-11.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-9.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-9.png>`
     - WGS
     -
     -
     - 0.11.0
     -
   * -
     - Tandem4
     - Indeterminate
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2D6-15.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-16.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-16.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA19719
     - 0.14.0
     -
   * - \*5, \*68+\*4
     - WholeDel1+Tandem1A
     - \*5/\*68+\*4
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2D6-13.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-11.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-11.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - HG01190
     - 0.4.1
     -
   * - \*2x2, \*68+\*4
     - WholeDup1+Tandem1A
     - \*2x2/\*68+\*4
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2D6-12.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-10.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-10.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA21781
     - 0.4.1
     -
   * -
     - ParalogPartialDel1
     - \*2/\*41
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2D6-14.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-15.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-15.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA19316
     - 0.13.0
     -
   * -
     - WholeDel1+Tandem3
     - Indeterminate
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2D6-17.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-20.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-20.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - HG03803
     - 0.16.0
     -
   * -
     - Unknown1
     - Indeterminate
     -
     -
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-14.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-14.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA18555
     - 0.12.0
     -
   * -
     - Unknown2
     - Indeterminate
     -
     -
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-18.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-18.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA19982
     - 0.14.0
     -

PyPGx was recently applied to the entire high-coverage WGS dataset from 1KGP
(N=2,504). Click `here <https://github.com/sbslee/1kgp-pgx-paper/blob/main/
sv-tables/CYP2D6.md>`__ to see individual SV calls for CYP2D6, and
corresponding copy number profiles and allele fraction profiles.

Phenotype summary for CYP2D6
----------------------------

Activity score is used for phenotype prediction.

.. list-table::
   :header-rows: 1

   * - Phenotype
     - Activity Score
     - Example
     - Priority
   * - Ultrarapid Metabolizer
     - 2.5 <= score
     - \*1/\*2x2
     - Abnormal/Priority/High Risk
   * - Normal Metabolizer
     - 1.25 <= score < 2.5
     - \*1/\*1
     - Normal/Routine/Low Risk
   * - Intermediate Metabolizer
     - 0.25 <= score < 1.25
     - \*1/\*4
     - Abnormal/Priority/High Risk
   * - Poor Metabolizer
     - 0 <= score < 0.25
     - \*4/\*5
     - Abnormal/Priority/High Risk
   * - Indeterminate
     - Unknown score
     - \*1/\*22
     - None

Recommendations for CYP2D6
--------------------------

- Amitriptyline

  "The CPIC Dosing Guideline update for amitriptyline recommends an
  alternative drug for CYP2D6 ultrarapid or poor metabolizers and CYP2C19
  ultrarapid, rapid or poor metabolizers. If amitriptyline is warranted,
  consider a 50% dose reduction in CYP2D6 or CYP2C19 poor metabolizers. For
  CYP2D6 intermediate metabolizers, a 25% dose reduction should be
  considered." (Source: `PharmGKB <https://www.pharmgkb.org
  /guidelineAnnotation/PA166105006>`__)

- Atomoxetine

  "The CPIC Dosing Guideline for atomoxetine provides therapeutic
  recommendations for CYP2D6 ultrarapid, normal, intermediate, and poor
  metabolizer, which includes guidance for plasma drug concentration testing,
  as a means to estimate atomoxetine exposure, if no clinical response and in
  the absence of adverse events after 2 weeks of therapy." (Source:
  `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/PA166181885>`__)

- Clomipramine

  "Tricyclic antidepressants have comparable pharmacokinetic properties, it
  may be reasonable to apply the CPIC Dosing Guideline for amitriptyline and
  CYP2C19, CYP2D6 to other tricyclics including clomipramine. The CPIC Dosing
  Guideline update for amitriptyline recommends an alternative drug for
  CYP2D6 ultrarapid or poor metabolizers and CYP2C19 ultrarapid, rapid or
  poor metabolizers. If amitriptyline is warranted, consider a 50% dose
  reduction in CYP2D6 or CYP2C19 poor metabolizers. For CYP2D6 intermediate
  metabolizers, a 25% dose reduction should be considered." (Source:
  `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/PA166105007>`__)

- Codeine

  "Alternate non-tramadol analgesics are recommended for CYP2D6 ultrarapid
  and poor metabolizers. A label recommended age- or weight-specific dose of
  codeine is warranted for CYP2D6 normal and intermediate metabolizers."
  (Source: `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/
  PA166104996>`__)

- Desipramine

  "Tricyclic antidepressants have comparable pharmacokinetic properties, it
  may be reasonable to apply the CPIC Dosing Guideline for
  amitriptyline/nortriptyline and CYP2C19, CYP2D6 to other tricyclics
  including
  desipramine. The CPIC Dosing Guideline update for nortriptyline recommends
  a 25% dose reduction for CYP2D6 intermediate metabolizers. For CYP2D6
  ultrarapid or poor metabolizers, an alternative drug should be considered.
  If nortriptyline is warranted, consider a 50% dose reduction in CYP2D6 poor
  metabolizers." (Source: `PharmGKB <https://www.pharmgkb.org/
  guidelineAnnotation/PA166105002>`__)

- Doxepin

  "Tricyclic antidepressants have comparable pharmacokinetic properties, it
  may be reasonable to apply the CPIC Dosing Guideline for amitriptyline and
  CYP2C19, CYP2D6 to other tricyclics including doxepin. The CPIC Dosing
  Guideline update for amitriptyline recommends an alternative drug for
  CYP2D6 ultrarapid or poor metabolizers and CYP2C19 ultrarapid, rapid or
  poor metabolizers. If amitriptyline is warranted, consider a 50% dose
  reduction in CYP2D6 or CYP2C19 poor metabolizers. For CYP2D6 intermediate
  metabolizers, a 25% dose reduction should be considered." (Source:
  `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/PA166105000>`__)

- Fluvoxamine

  "The CPIC Dosing Guideline for the selective serotonin reuptake inhibitor
  fluvoxamine recommends to consider a 25-50% reduction of recommended
  starting dose and titrate to response or use an alternative drug not
  metabolized by CYP2D6 for CYP2D6 poor metabolizers." (Source:
  `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/PA166127637>`__)

- Hydrocodone

  "CYP2D6 intermediate and poor metabolizers should initiate hydrocodone
  therapy using the label recommended age- or weight-specific dosing.
  However, if there is no response to hydrocodone in these patients, an
  alternative analgesic should be considered.

  There is insufficient evidence to provide a dosing recommendation for
  hydrocodone in CYP2D6 ultrarapid metabolizers." (Source: `PharmGKB
  <https://www.pharmgkb.org/guidelineAnnotation/PA166228121>`__)

- Imipramine

  "Tricyclic antidepressants have comparable pharmacokinetic properties, it
  may be reasonable to apply the CPIC Dosing Guideline for amitriptyline and
  CYP2C19, CYP2D6 to other tricyclics including imipramine. The CPIC Dosing
  Guideline update for amitriptyline recommends an alternative drug for
  CYP2D6 ultrarapid or poor metabolizers and CYP2C19 ultrarapid, rapid or
  poor metabolizers. If amitriptyline is warranted, consider a 50% dose
  reduction in CYP2D6 or CYP2C19 poor metabolizers. For CYP2D6 intermediate
  metabolizers, a 25% dose reduction should be considered." (Source:
  `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/PA166104999>`__)

- Nortriptyline

  "The CPIC Dosing Guideline update for nortriptyline recommends a 25% dose
  reduction for CYP2D6 intermediate metabolizers. For CYP2D6 ultrarapid or
  poor metabolizers, an alternative drug should be considered. If
  nortriptyline is warranted, consider a 50% dose reduction in CYP2D6 poor
  metabolizers." (Source: `PharmGKB <https://www.pharmgkb.org/
  guidelineAnnotation/PA166104999>`__)

- Ondansetron

  "The CPIC dosing guideline for ondansetron recommends selecting an
  alternate drug for CYP2D6 ultrarapid metabolizers. It is recommended that
  the alternate drug not be predominantly metabolized by CYP2D6 (eg.
  granisetron)." (Source: `PharmGKB <https://www.pharmgkb.org/
  guidelineAnnotation/PA166161954>`__)

- Paroxetine

  "The CPIC Dosing Guideline for the selective serotonin reuptake inhibitor
  paroxetine recommends an alternative drug not predominantly metabolized by
  CYP2D6 for CYP2D6 ultrarapid metabolizers and for CYP2D6 poor metabolizers.
  For CYP2D6 poor metabolizers, if paroxetine use is warranted, consider a
  50% reduction of recommended starting dose and titrate to response."
  (Source: `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/
  PA166127636>`__)

- Tamoxifen

  "The CPIC Dosing Guideline for tamoxifen recommends the use of alternative
  hormonal therapy such as an aromatase inhibitor for postmenopausal women or
  aromatase inhibitor along with ovarian function suppression in
  premenopausal women for CYP2D6 poor metabolizer, if aromatase inhibitor use
  is not contraindicated. For CYP2D6 intermediate metabolizers and CYP2D6
  allele combinations resulting in an activity score (AS) of 1 the
  recommendation is to consider the recommendations stated for the CYP2D6
  poor metabolizer. If aromatase inhibitor use is contraindicated,
  consideration should be given to use a higher but FDA approved tamoxifen
  dose for CYP2D6 intermediate metabolizers and CYP2D6 allele combinations
  resulting in an AS of 1. For poor metabolizer, higher dose tamoxifen (40
  mg/day) increases but does not normalize endoxifen concentrations and can
  be considered if there are contraindications to aromatase inhibitor
  therapy." (Source: `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation
  /PA166176068>`__)

- Tramadol

  "Alternate non-codeine analgesics are recommended for CYP2D6 ultrarapid and
  poor metabolizers. A label recommended age- or weight-specific dose of
  tramadol is warranted for CYP2D6 normal and intermediate metabolizers."
  (Source: `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/
  PA166228101>`__)

- Trimipramine

  "Tricyclic antidepressants have comparable pharmacokinetic properties, it
  may be reasonable to apply the CPIC Dosing Guideline for amitriptyline and
  CYP2C19, CYP2D6 to other tricyclics including trimipramine. The CPIC Dosing
  Guideline update for amitriptyline recommends an alternative drug for
  CYP2D6 ultrarapid or poor metabolizers and CYP2C19 ultrarapid, rapid or
  poor metabolizers. If amitriptyline is warranted, consider a 50% dose
  reduction in CYP2D6 or CYP2C19 poor metabolizers. For CYP2D6 intermediate
  metabolizers, a 25% dose reduction should be considered." (Source:
  `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/PA166105001>`__)

- Tropisetron

  "The CPIC dosing guideline for tropisetron recommends selecting an
  alternate drug for CYP2D6 ultrarapid metabolizers. It is recommended that
  the alternate drug not be predominantly metabolized by CYP2D6 (eg.
  granisetron)." (Source: `PharmGKB <https://www.pharmgkb.org/
  guidelineAnnotation/PA166161955>`__)

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

- `PharmVar CYP2D6 page <https://www.pharmvar.org/gene/CYP2D6>`__
- `CPIC® Guideline for Tamoxifen based on CYP2D6 genotype <https://cpicpgx.org/guidelines/cpic-guideline-for-tamoxifen-based-on-cyp2d6-genotype/>`__

CYP2E1
======

SV summary for CYP2E1
---------------------

Below is comprehensive summary of SV described from real NGS studies:

.. list-table::
   :header-rows: 1

   * - SV Alleles
     - SV Name
     - Genotype
     - Reference
     - Gene Model
     - GRCh37
     - GRCh38
     - Data Type
     - Source
     - Coriell ID
     - Version
     - Description
   * -
     - Normal
     - \*1/\*7
     - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2E1-1.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2E1-5.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2E1-5.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA10831
     - 0.4.1
     -
   * -
     - WholeDel1
     - Indeterminate
     -
     -
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2E1-9.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2E1-9.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - HG03445
     - 0.16.0
     -
   * - \*1x2
     - WholeDup1
     - \*1/\*1x2
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2E1-3.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2E1-4.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2E1-4.png>`
     - WGS
     -
     -
     - 0.4.1
     -
   * - \*7x2
     - WholeDup1
     - \*1/\*7x2
     - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2E1-3.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2E1-2.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2E1-2.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA19095
     - 0.4.1
     -
   * - \*1x2
     - WholeDup2
     - \*1x2/\*7
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2E1-3.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2E1-6.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2E1-6.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA19225
     - 0.12.0
     -
   * - \*S1
     - PartialDup1
     - \*1/\*S1
     - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2E1-2.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2E1-1.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2E1-1.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA19920
     - 0.4.1
     - \*S1 is linked to \*7.
   * - \*S1
     - PartialDup1Hom
     - \*S1/\*S1
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2E1-5.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2E1-7.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2E1-7.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA19309
     - 0.13.0
     -
   * - \*7x3
     - WholeMultip1
     - \*7/\*7x3
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2E1-4.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2E1-3.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2E1-3.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA19908
     - 0.4.1
     -
   * -
     - WholeMultip2
     - Indeterminate
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2E1-6.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2E1-8.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2E1-8.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA20291
     - 0.14.0
     -
   * -
     - WholeDup1+PartialDup1
     - Indeterminate
     -
     -
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2E1-10.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2E1-10.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - HG03401
     - 0.16.0
     -

PyPGx was recently applied to the entire high-coverage WGS dataset from 1KGP
(N=2,504). Click `here <https://github.com/sbslee/1kgp-pgx-paper/blob/main/
sv-tables/CYP2E1.md>`__ to see individual SV calls for CYP2E1, and
corresponding copy number profiles and allele fraction profiles.

Resources for CYP2E1
--------------------

- `PharmVar CYP2E1 page <https://www.pharmvar.org/gene/CYP2E1>`__

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
     - Priority
   * - Normal Metabolizer
     - \*1/\*1
     - Abnormal/Priority/High Risk
   * - Intermediate Metabolizer
     - \*1/\*3
     - Abnormal/Priority/High Risk
   * - Possible Intermediate Metabolizer
     - \*1/\*2
     - Abnormal/Priority/High Risk
   * - Poor Metabolizer
     - \*6/\*6
     - Normal/Routine/Low Risk
   * - Indeterminate
     - \*2/\*2
     - None

Recommendations for CYP3A5
--------------------------

- Tacrolimus

  "The CPIC dosing guideline for tacrolimus recommends increasing the starting
  dose by 1.5 to 2 times the recommended starting dose in patients who are
  CYP3A5 intermediate or extensive metabolizers, though total starting dose
  should not exceed 0.3 mg/kg/day. Therapeutic drug monitoring should also
  be used to guide dose adjustments." (Source: `PharmGKB <https://
  www.pharmgkb.org/guidelineAnnotation/PA166124619>`__)

Resources for CYP3A5
--------------------

- `PharmVar CYP3A5 page <https://www.pharmvar.org/gene/CYP3A5>`__
- `CPIC® Guideline for Tacrolimus and CYP3A5 <https://cpicpgx.org/guidelines/guideline-for-tacrolimus-and-cyp3a5/>`__
- `PharmGKB: Annotation of CPIC Guideline for tacrolimus and CYP3A5 <https://www.pharmgkb.org/guidelineAnnotation/PA166124619>`__

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

  * - SV Alleles
    - SV Name
    - Genotype
    - Reference
    - Gene Model
    - GRCh37
    - GRCh38
    - Data Type
    - Source
    - Coriell ID
    - Version
    - Description
  * -
    - Normal
    - \*1/\*3
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP4F2-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP4F2-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP4F2-2.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - HG01190
    - 0.11.0
    -
  * - \*DEL
    - WholeDel1
    - \*1/\*DEL
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP4F2-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP4F2-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP4F2-1.png>`
    - WGS
    -
    -
    - 0.11.0
    -

PyPGx was recently applied to the entire high-coverage WGS dataset from 1KGP
(N=2,504). Click `here <https://github.com/sbslee/1kgp-pgx-paper/blob/main/
sv-tables/CYP4F2.md>`__ to see individual SV calls for CYP4F2, and
corresponding copy number profiles and allele fraction profiles.

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
     - Priority
   * - Normal Metabolizer
     - 2 == score
     - Reference/Reference
     - Normal/Routine/Low Risk
   * - Intermediate Metabolizer
     - 1 <= score < 2
     - Reference/c.1905+1G>A (\*2A)
     - Abnormal/Priority/High Risk
   * - Poor Metabolizer
     - 0 <= score < 1
     - c.295_298delTCAT (\*7)/c.703C>T (\*8)
     - Abnormal/Priority/High Risk

Recommendations for DPYD
------------------------

- Capecitabine

  "The CPIC Dosing Guideline for 5-fluorouracil and capecitabine recommends
  an alternative drug for patients who are DPYD poor metabolizers with an
  activity score of 0. In those who are poor metabolizers with an activity
  score of 0.5, an alternative drug is also recommended, but if this is not
  considered a suitable therapeutic option, 5-fluorouracil or capecitabine
  should be administered at a strongly reduced dose with early therapeutic
  drug monitoring. Patients who are intermediate metabolizers with an
  activity score of 1 or 1.5 should receive a dose reduction of 50%. Patients
  with the c.[2846A>T];[2846A>T] genotype may require a >50% dose reduction."
  (Source: `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/
  PA166109594>`__)

- Fluorouracil

  "The CPIC Dosing Guideline for 5-fluorouracil and capecitabine recommends
  an alternative drug for patients who are DPYD poor metabolizers with an
  activity score of 0. In those who are poor metabolizers with an activity
  score of 0.5, an alternative drug is also recommended, but if this is not
  considered a suitable therapeutic option, 5-fluorouracil or capecitabine
  should be administered at a strongly reduced dose with early therapeutic
  drug monitoring. Patients who are intermediate metabolizers with an
  activity score of 1 or 1.5 should receive a dose reduction of 50%. Patients
  with the c.[2846A>T];[2846A>T] genotype may require a >50% dose reduction."
  (Source: `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/
  PA166122686>`__)

Resources for DPYD
------------------

- `PharmVar DPYD page <https://www.pharmvar.org/gene/DPYD>`__
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
      - Priority
    * - Favorable Response
      - Reference/Reference
      - None
    * - Unfavorable Response
      - Reference/Leiden
      - None

Recommendations for F5
----------------------

- Hormonal contraceptives

  "In individuals who carry the Factor V Leiden allele (rs6025 T) and have a
  family history of thrombotic events, estrogen-containing oral
  contraceptives should be avoided and alternative forms of contraception
  used." (Source: `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/
  PA166104955>`__)

Resources for F5
----------------

- `PharmGKB: Annotation of DPWG Guideline for hormonal contraceptives for systemic use and F5 <https://www.pharmgkb.org/chemical/PA452637/guidelineAnnotation/PA166104955>`__

G6PD
====

SV summary for G6PD
-------------------

Since the gene is located on X chromosome, its copy number differs between
females (N=2) and males (N=1). Technically speaking, this difference is not a
SV event, but it is treated as such by PyPGx for genotyping purposes (i.e.
sex determination).

Below is comprehensive summary of SV described from real NGS studies:

.. list-table::
  :header-rows: 1

  * - SV Alleles
    - SV Name
    - Genotype
    - Reference
    - Gene Model
    - GRCh37
    - GRCh38
    - Data Type
    - Source
    - Coriell ID
    - Version
    - Description
  * -
    - Female
    - B (reference)/B (reference)
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-G6PD-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-G6PD-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-G6PD-1.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - HG00276
    - 0.12.0
    -
  * - MALE
    - Male
    - MALE/B (reference)
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-G6PD-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-G6PD-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-G6PD-2.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - HG00436
    - 0.12.0
    -

PyPGx was recently applied to the entire high-coverage WGS dataset from 1KGP
(N=2,504). Click `here <https://github.com/sbslee/1kgp-pgx-paper/blob/main/
sv-tables/G6PD.md>`__ to see individual SV calls for G6PD, and
corresponding copy number profiles and allele fraction profiles.

Phenotype summary for G6PD
--------------------------

Diplotype-phenotype mapping is used for phenotype prediction.

 .. list-table::
    :header-rows: 1

    * - Phenotype
      - Example
      - Priority
    * - G6PD Normal
      - A/B (reference)
      - Normal Risk
    * - G6PD Deficient with CNSHA
      - Alhambra/Bari
      - Abnormal/Priority/High Risk
    * - G6PD Variable
      - 202G>A_376A>G_1264C>G/B (reference)
      - Abnormal/Priority/High Risk
    * - G6PD Deficient
      - 202G>A_376A>G_1264C>G/A- 680T_376G
      - Abnormal/Priority/High Risk
    * - G6PD Indeterminate
      - A- 202A_376G/Dagua
      - Abnormal/Priority/High Risk

GSTM1
=====

SV summary for GSTM1
--------------------

This gene is known to have an extremely high rate of gene deletion
polymorphism in the population and thus requires SV analysis.

Below is comprehensive summary of SV described from real NGS studies:

.. list-table::
  :header-rows: 1

  * - SV Alleles
    - SV Name
    - Genotype
    - Reference
    - Gene Model
    - GRCh37
    - GRCh38
    - Data Type
    - Source
    - Coriell ID
    - Version
    - Description
  * -
    - Normal
    - \*A/\*B
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-GSTM1-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTM1-5.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTM1-5.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA06991
    - 0.4.1
    -
  * - \*0
    - WholeDel1
    - \*0/\*A
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-GSTM1-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTM1-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTM1-1.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA18855
    - 0.4.1
    -
  * - \*0
    - WholeDel1Hom
    - \*0/\*0
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-GSTM1-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTM1-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTM1-2.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA10831
    - 0.4.1
    -
  * - \*0
    - WholeDel2
    - \*0/\*A
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-GSTM1-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTM1-10.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTM1-10.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA21097
    - 0.15.0
    -
  * - \*Ax2
    - WholeDup1
    - \*A/\*Ax2
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-GSTM1-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTM1-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTM1-3.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA19908
    - 0.4.1
    -
  * - \*Bx2
    - WholeDup1
    - \*A/\*Bx2
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-GSTM1-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTM1-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTM1-4.png>`
    - WGS
    -
    -
    - 0.4.1
    -
  * -
    - NoncodingDel1
    - \*A/\*B
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-GSTM1-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTM1-6.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTM1-6.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - NA19005
    - 0.12.0
    -
  * - \*0
    - WholeDel1+NoncodingDel1
    - \*0/\*A
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-GSTM1-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTM1-7.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTM1-7.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - NA06984
    - 0.12.0
    -
  * -
    - PartialDup1
    - Indeterminate
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-GSTM1-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTM1-8.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTM1-8.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - NA19908
    - 0.14.0
    -
  * -
    - WholeDel1+WholeDel2
    - Indeterminate
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-GSTM1-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTM1-9.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTM1-9.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - NA20506
    - 0.14.0
    -

PyPGx was recently applied to the entire high-coverage WGS dataset from 1KGP
(N=2,504). Click `here <https://github.com/sbslee/1kgp-pgx-paper/blob/main/
sv-tables/GSTM1.md>`__ to see individual SV calls for GSTM1, and
corresponding copy number profiles and allele fraction profiles.

GSTT1
=====

GRCh38 data for GSTT1
---------------------

GSTT1 is located on ``chr22`` for GRCh37 but on ``chr22_KI270879v1_alt``
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

  * - SV Alleles
    - SV Name
    - Genotype
    - Reference
    - Gene Model
    - GRCh37
    - GRCh38
    - Data Type
    - Source
    - Coriell ID
    - Version
    - Description
  * -
    - Normal
    - \*A/\*A
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-GSTT1-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTT1-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTT1-3.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA07055
    - 0.4.1
    -
  * - \*0
    - WholeDel1
    - \*0/\*A
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-GSTT1-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTT1-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTT1-1.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA19908
    - 0.4.1
    -
  * - \*0
    - WholeDel1Hom
    - \*0/\*0
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-GSTT1-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTT1-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTT1-2.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA11832
    - 0.4.1
    -

IFNL3
=====

Phenotype summary for IFNL3
---------------------------

Diplotype-phenotype mapping is used for phenotype prediction.

.. list-table::
   :header-rows: 1

   * - Phenotype
     - Example
     - Priority
   * - Favorable Response
     - Reference/Reference
     - None
   * - Unfavorable Response
     - Reference/rs12979860
     - None
   * - Indeterminate
     - Reference/rs8099917
     - None

Recommendations for IFNL3
-------------------------

- Peginterferon alfa-2a

"IFNL3 (IL28B) variation (rs12979860) is the strongest baseline predictor of
response to PEG-interferon-alpha-containing regimens in HCV genotype 1
patients. Patients with the favorable response genotype (rs12979860 CC) have
increased likelihood of response (higher SVR rate) to PEG-interferon-alpha-
containing regimens as compared to patients with unfavorable response
genotype (rs12979860 CT or TT). Consider implications before initiating
PEG-IFN alpha and RBV containing regimens." (Source: `PharmGKB <https://www.
pharmgkb.org/chemical/PA164749390/guidelineAnnotation/PA166110235>`__)

- Peginterferon alfa-2b

"IFNL3 (IL28B) variation (rs12979860) is the strongest baseline predictor of
response to PEG-interferon-alpha-containing regimens in HCV genotype 1
patients. Patients with the favorable response genotype (rs12979860 CC) have
increased likelihood of response (higher SVR rate) to PEG-interferon-alpha-
containing regimens as compared to patients with unfavorable response
genotype (rs12979860 CT or TT). Consider implications before initiating
PEG-IFN alpha and RBV containing regimens." (Source: `PharmGKB <https://www.
pharmgkb.org/chemical/PA164784024/guidelineAnnotation/PA166110235>`__)

- Ribavirin

"IFNL3 (IL28B) variation (rs12979860) is the strongest baseline predictor of
response to PEG-interferon-alpha-containing regimens in HCV genotype 1
patients. Patients with the favorable response genotype (rs12979860 CC) have
increased likelihood of response (higher SVR rate) to PEG-interferon-alpha-
containing regimens as compared to patients with unfavorable response
genotype (rs12979860 CT or TT). Consider implications before initiating
PEG-IFN alpha and RBV containing regimens." (Source: `PharmGKB <https://www.
pharmgkb.org/chemical/PA451241/guidelineAnnotation/PA166110235>`__)

Resources for IFNL3
-------------------

- `PharmGKB: Annotation of CPIC Guideline for peginterferon alfa-2a,peginterferon alfa-2b,ribavirin and IFNL3 <https://www.pharmgkb.org/guidelineAnnotation/PA166110235>`__
- `CPIC® Guideline for PEG Interferon-Alpha-Based Regimens and IFNL3 <https://cpicpgx.org/guidelines/guideline-for-peg-interferon-alpha-based-regimens-and-ifnl3/>`__

MT-RNR1
=======

Phenotype summary for MT-RNR1
-----------------------------

Diplotype-phenotype mapping is used for phenotype prediction.

.. list-table::
   :header-rows: 1

   * - Phenotype
     - Example
     - Priority
   * - Normal Risk of Aminoglycoside-Induced Hearing Loss
     - Reference/Reference
     - None
   * - Increased Risk of Aminoglycoside-Induced Hearing Loss
     - m.1494C>T/m.1494C>T
     - None
   * - Uncertain Risk of Aminoglycoside-Induced Hearing Loss
     - m.663A>G/m.663A>G
     - None

NUDT15
======

Phenotype summary for NUDT15
----------------------------

Diplotype-phenotype mapping is used for phenotype prediction.

.. list-table::
   :header-rows: 1

   * - Phenotype
     - Example
     - Priority
   * - Normal Metabolizer
     - \*1/\*1
     - Normal/Routine/Low risk
   * - Intermediate Metabolizer
     - \*1/\*2
     - Abnormal/Priority/High Risk
   * - Possible Intermediate Metabolizer
     - \*3/\*4
     - Abnormal/Priority/High Risk
   * - Poor Metabolizer
     - \*2/\*3
     - Abnormal/Priority/High Risk
   * - Indeterminate
     - \*1/\*4
     - Abnormal/Priority/High Risk

Recommendations for NUDT15
--------------------------

- Azathioprine

  "Consider an alternate agent or extreme dose reduction of azathioprine for
  patients who are TPMT or NUDT15 poor metabolizers. Start at 30-80% of
  target dose for patients who are TPMT or NUDT15 intermediate metabolizers."
  (Source: `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/
  PA166104933>`__)

- Mercaptopurine

  "Consider an alternate agent or extreme dose reduction of mercaptopurine
  for patients who are TPMT or NUDT15 poor metabolizers. Start at 30-80% of
  target dose for patients who are TPMT or NUDT15 intermediate metabolizers."
  (Source: `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/
  PA166104945>`__)

- Thioguanine

  "Consider an alternate agent or extreme dose reduction of thioguanine for
  patients who are TPMT or NUDT15 poor metabolizers. Start at 50-80% of
  target dose for patients who are TPMT or NUDT15 intermediate metabolizers."
  (Source: `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/
  PA166104965>`__)

Resources for NUDT15
--------------------

- `PharmVar NUDT15 page <https://www.pharmvar.org/gene/NUDT15>`__
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
    - Priority
  * - Uncertain Susceptibility
    - Reference/Reference
    - Normal Risk
  * - Malignant Hyperthermia Susceptibility
    - Reference/c.103T>C
    - Abnormal/Priority/High Risk

Recommendations for RYR1
------------------------

- Desflurane

  "The CPIC Dosing Guideline recommends that halogenated volatile anesthetics
  uch as desflurane, enflurane, halothane, isoflurane, methoxyflurane,
  sevoflurane and the depolarizing muscle relaxants succinylcholine are
  relatively contraindicated in persons with malignant hyperthermia
  susceptibility (MHS). See full guideline for disclaimers, further details
  and supporting evidence." (Source: `PharmGKB <https://www.pharmgkb.org/
  chemical/PA164749136/guidelineAnnotation/PA166180457>`__)

- Enflurane

  "The CPIC Dosing Guideline recommends that halogenated volatile anesthetics
  such as desflurane, enflurane, halothane, isoflurane, methoxyflurane,
  sevoflurane and the depolarizing muscle relaxants succinylcholine are
  relatively contraindicated in persons with malignant hyperthermia
  susceptibility (MHS). See full guideline for disclaimers, further details
  and supporting evidence." (Source: `PharmGKB <https://www.pharmgkb.org/
  chemical/PA449461/guidelineAnnotation/PA166180457>`__)

- Halothane

  "The CPIC Dosing Guideline recommends that halogenated volatile anesthetics
  such as desflurane, enflurane, halothane, isoflurane, methoxyflurane,
  sevoflurane and the depolarizing muscle relaxants succinylcholine are
  relatively contraindicated in persons with malignant hyperthermia
  susceptibility (MHS). See full guideline for disclaimers, further details
  and supporting evidence." (Source: `PharmGKB <https://www.pharmgkb.org/
  chemical/PA449845/guidelineAnnotation/PA166180457>`__)

- Isoflurane

  "The CPIC Dosing Guideline recommends that halogenated volatile anesthetics
  such as desflurane, enflurane, halothane, isoflurane, methoxyflurane,
  sevoflurane and the depolarizing muscle relaxants succinylcholine are
  relatively contraindicated in persons with malignant hyperthermia
  susceptibility (MHS). See full guideline for disclaimers, further details
  and supporting evidence." (Source: `PharmGKB <https://www.pharmgkb.org/
  chemical/PA450106/guidelineAnnotation/PA166180457>`__)

- Methoxyflurane

  "The CPIC Dosing Guideline recommends that halogenated volatile anesthetics
  such as desflurane, enflurane, halothane, isoflurane, methoxyflurane,
  sevoflurane and the depolarizing muscle relaxants succinylcholine are
  relatively contraindicated in persons with malignant hyperthermia
  susceptibility (MHS). See full guideline for disclaimers, further details
  and supporting evidence." (Source: `PharmGKB <https://www.pharmgkb.org/
  chemical/PA450434/guidelineAnnotation/PA166180457>`__)

- Sevoflurane

  "The CPIC Dosing Guideline recommends that halogenated volatile anesthetics
  such as desflurane, enflurane, halothane, isoflurane, methoxyflurane,
  sevoflurane and the depolarizing muscle relaxants succinylcholine are
  relatively contraindicated in persons with malignant hyperthermia
  susceptibility (MHS). See full guideline for disclaimers, further details
  and supporting evidence." (Source: `PharmGKB <https://www.pharmgkb.org/
  chemical/PA451341/guidelineAnnotation/PA166180457>`__)

- Succinylcholine

  "The CPIC Dosing Guideline recommends that halogenated volatile anesthetics
  such as desflurane, enflurane, halothane, isoflurane, methoxyflurane,
  sevoflurane and the depolarizing muscle relaxants succinylcholine are
  relatively contraindicated in persons with malignant hyperthermia
  susceptibility (MHS). See full guideline for disclaimers, further details
  and supporting evidence." (Source: `PharmGKB <https://www.pharmgkb.org/
  chemical/PA451522/guidelineAnnotation/PA166180457>`__)

Resources for RYR1
------------------

- `PharmGKB: Annotation of CPIC Guideline for desflurane and CACNA1S, RYR1 <https://www.pharmgkb.org/chemical/PA164749136/guidelineAnnotation/PA166180457>`__
- `CPIC® Guideline for Potent Volatile Anesthetic Agents and Succinylcholine and RYR1 and CACNA1S <https://cpicpgx.org/guidelines/cpic-guideline-for-ryr1-and-cacna1s/>`__
- `Clinical Pharmacogenetics Implementation Consortium (CPIC) Guideline for the Use of Potent Volatile Anesthetic Agents and Succinylcholine in the Context of RYR1 or CACNA1S Genotypes <https://doi.org/10.1002/cpt.1319>`__

SLC22A2
=======

SV summary for SLC22A2
----------------------

Below is comprehensive summary of SV described from real NGS studies:

.. list-table::
  :header-rows: 1

  * - SV Alleles
    - SV Name
    - Genotype
    - Reference
    - Gene Model
    - GRCh37
    - GRCh38
    - Data Type
    - Source
    - Coriell ID
    - Version
    - Description
  * -
    - Normal
    - \*1/\*3
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-SLC22A2-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SLC22A2-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SLC22A2-1.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - HG01190
    - 0.4.1
    -
  * - \*S1
    - NoncodingDel1
    - \*1/\*S1
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-SLC22A2-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SLC22A2-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SLC22A2-2.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA18855
    - 0.4.1
    -
  * - \*S1
    - NoncodingDel1Hom
    - \*S1/\*S1
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-SLC22A2-6.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SLC22A2-6.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SLC22A2-6.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - HG02337
    - 0.16.0
    -
  * - \*S2
    - PartialDel1
    - \*1/\*S2
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-SLC22A2-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SLC22A2-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SLC22A2-3.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA19819
    - 0.4.1
    -
  * - \*S1, \*S2
    - NoncodingDel1+PartialDel1
    - \*S1/\*S2
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-SLC22A2-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SLC22A2-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SLC22A2-4.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - NA19030
    - 0.13.0
    -
  * -
    - PartialDup1
    - Indeterminate
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-SLC22A2-5.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SLC22A2-5.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SLC22A2-5.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - NA20813
    - 0.14.0
    -

PyPGx was recently applied to the entire high-coverage WGS dataset from 1KGP
(N=2,504). Click `here <https://github.com/sbslee/1kgp-pgx-paper/blob/main/
sv-tables/SLC22A2.md>`__ to see individual SV calls for SLC22A2, and
corresponding copy number profiles and allele fraction profiles.

SLCO1B1
=======

Phenotype summary for SLCO1B1
-----------------------------

Diplotype-phenotype mapping is used for phenotype prediction.

.. list-table::
   :header-rows: 1

   * - Phenotype
     - Example
     - Priority
   * - Increased Function
     - \*14/\*14
     - None
   * - Normal Function
     - \*1/\*1
     - Normal/Routine/Low Risk
   * - Possible Decreased Function
     - \*2/\*15
     - Abnormal/Priority/High Risk
   * - Decreased Function
     - \*1/\*5
     - Abnormal/Priority/High Risk
   * - Poor Function
     - \*5/\*5
     - Abnormal/Priority/High Risk
   * - Indeterminate
     - \*2/\*38
     - None

Recommendations for SLCO1B1
---------------------------

- Atorvastatin

  "Prescribe ≤20mg for patients with SLCO1B1 poor function phenotype and
  ≤40mg for patients with SLCO1B1 decreased or possible decreased phenotype
  as a starting dose. Adjust doses of atorvastatin based on disease-specific
  guidelines. Prescriber should be aware of possible increased risk for
  myopathy especially for 40mg dose." (Source: `PharmGKB <https://www.
  pharmgkb.org/guidelineAnnotation/PA166262221>`__)

- Fluvastatin

  "CYP2C9 IMs should avoid fluvastatin doses greater than 40mg while CYP2C9
  PMs should avoid doses greater than 20mg. If higher doses are required for
  desired efficacy, an alternative statin should be considered. Patients with
  SLCO1B1 poor function should also avoid fluvastatin doses greater than 40mg
  and and adjust doses of fluvastatin based on disease-specific guidelines.
  Patients with both SLCO1B1 poor function and CYP2C9 IM/PM should be
  prescribed an alternative statin depending on the desired potency."
  (Source: `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/
  PA166262341>`__)

- Lovastatin

  "Prescribe an alternative statin depending on the desired potency for
  patients with SLCO1B1 decreased function, possible decreased function or
  poor function phenotype. If lovastatin therapy is warranted in patients
  with SLCO1B1 decreased or possible decreased phenotype, limit dose to
  <20mg/day." (Source: `PharmGKB <https://www.pharmgkb.org/
  guidelineAnnotation/PA166262241>`__)

- Pitavastatin

  "Prescribe ≤1mg as a starting dose for patients with SLCO1B1 poor function
  phenotype. Prescribe ≤2mg as a starting dose for patients with SLCO1B1
  decreased or possible decreased phenotype. Adjust doses of pitavastatin
  based on disease-specific guidelines. Consider an alternative statin or
  combination therapy if higher doses are needed." (Source: `PharmGKB
  <https://www.pharmgkb.org/guidelineAnnotation/PA166262261>`__)

- Pravastatin

  "Prescribe ≤40mg as a starting dose and adjust doses of pravastatin based
  on disease-specific guidelines for patients with SLCO1B1 poor function
  phenotype. Prescribe desired starting dose and adjust doses of pravastatin
  based on disease-specific guidelines for patients with SLCO1B1 decreased or
  possible decreased phenotype. Prescriber should be aware of possible
  increased risk for myopathy with pravastatin especially with doses >40mg
  per day." (Source: `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation
  /PA166262281>`__)

- Rosuvastatin

  "Prescribe ≤20mg as a starting dose and adjust doses of rosuvastatin based
  on disease-specific and specific population guidelines for patients who are
  SLCO1B1 or ABCG2 poor function phenotype. If dose >20mg needed for desired
  efficacy, consider combination therapy (i.e. rosuvastatin plus non-statin
  guideline directed medical therapy). Patients with both ABCG2 poor function
  and SLCO1B1 poor/decreased function should be prescribed ≤10mg as a
  starting dose." (Source: `PharmGKB <https://www.pharmgkb.org/
  guidelineAnnotation/PA166262321>`__)

- Simvastatin

  "Prescribe an alternative statin depending on the desired potency for
  patients with SLCO1B1 decreased function, possible decreased function or
  poor function phenotype. If simvastatin therapy is warranted in patients
  with SLCO1B1 decreased or possible decreased phenotype, limit dose to
  <20mg/day." (Source: `PharmGKB <https://www.pharmgkb.org/
  guidelineAnnotation/PA166105005>`__)

Resources for SLCO1B1
---------------------

- `PharmVar SLCO1B1 page <https://www.pharmvar.org/gene/SLCO1B1>`__
- `CPIC® Guideline for Simvastatin and SLCO1B1 <https://cpicpgx.org/guidelines/guideline-for-simvastatin-and-slco1b1/>`__
- `The Clinical Pharmacogenetics Implementation Consortium Guideline for SLCO1B1, ABCG2, and CYP2C9 genotypes and Statin-Associated Musculoskeletal Symptoms <https://ascpt.onlinelibrary.wiley.com/doi/10.1002/cpt.2557>`__
- `PharmGKB: Annotation of CPIC Guideline for rosuvastatin and ABCG2, SLCO1B1 <https://www.pharmgkb.org/guidelineAnnotation/PA166262321>`__

SULT1A1
=======

SV summary for SULT1A1
----------------------

Below is comprehensive summary of SV described from real NGS studies:

.. list-table::
  :header-rows: 1

  * - SV Alleles
    - SV Name
    - Genotype
    - Reference
    - Gene Model
    - GRCh37
    - GRCh38
    - Data Type
    - Source
    - Coriell ID
    - Version
    - Description
  * -
    - Normal
    - \*1/\*2
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-SULT1A1-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SULT1A1-6.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SULT1A1-6.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA06991
    - 0.11.0
    -
  * - \*DEL
    - WholeDel1
    - \*1/\*DEL
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-SULT1A1-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SULT1A1-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SULT1A1-1.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA18942
    - 0.11.0
    -
  * - \*DEL
    - WholeDel1Hom
    - \*DEL/\*DEL
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-SULT1A1-7.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SULT1A1-7.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SULT1A1-7.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - NA20874
    - 0.14.0
    -
  * - \*1x2
    - WholeDup1
    - \*1x2/\*2
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-SULT1A1-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SULT1A1-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SULT1A1-2.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA18509
    - 0.11.0
    -
  * - \*1x3
    - WholeMultip1
    - \*1x3/\*2
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-SULT1A1-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SULT1A1-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SULT1A1-3.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA18868
    - 0.11.0
    -
  * - \*1x4
    - WholeMultip2
    - \*1x4/\*2
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-SULT1A1-5.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SULT1A1-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SULT1A1-4.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA18484
    - 0.11.0
    -
  * - \*1x3, \*2x2
    - WholeMultip2
    - \*1x3/\*2x2
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-SULT1A1-6.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SULT1A1-5.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SULT1A1-5.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA19143
    - 0.11.0
    -
  * -
    - Unknown1
    - Indeterminate
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-SULT1A1-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SULT1A1-8.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SULT1A1-8.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - HG01085
    - 0.15.0
    -
  * -
    - Unknown2
    - Indeterminate
    -
    -
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SULT1A1-9.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SULT1A1-9.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    -
    - 0.16.0
    -
  * -
    - Unknown3
    - Indeterminate
    -
    -
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SULT1A1-10.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SULT1A1-10.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - HG03854
    - 0.16.0
    -
  * -
    - Unknown4
    - Indeterminate
    -
    -
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SULT1A1-11.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SULT1A1-11.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - HG03742
    - 0.16.0
    -

PyPGx was recently applied to the entire high-coverage WGS dataset from 1KGP
(N=2,504). Click `here <https://github.com/sbslee/1kgp-pgx-paper/blob/main/
sv-tables/SULT1A1.md>`__ to see individual SV calls for SULT1A1, and
corresponding copy number profiles and allele fraction profiles.

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
     - Priority
   * - Normal Metabolizer
     - \*1/\*1
     - Normal/Routine/Low Risk
   * - Possible Intermediate Metabolizer
     - \*3A/\*12
     - Abnormal/Priority/High Risk
   * - Intermediate Metabolizer
     - \*1/\*2
     - Abnormal/Priority/High Risk
   * - Poor Metabolizer
     - \*2/\*3A
     - Abnormal/Priority/High Risk
   * - Indeterminate
     - \*1/\*18
     - Abnormal/Priority/High Risk

Recommendations for TPMT
------------------------

- Azathioprine

  "Consider an alternate agent or extreme dose reduction of azathioprine for
  patients who are TPMT or NUDT15 poor metabolizers. Start at 30-80% of
  target dose for patients who are TPMT or NUDT15 intermediate metabolizers."
  (Source: `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/
  PA166104933>`__)

- Mercaptopurine

  "Consider an alternate agent or extreme dose reduction of mercaptopurine
  for patients who are TPMT or NUDT15 poor metabolizers. Start at 30-80% of
  target dose for patients who are TPMT or NUDT15 intermediate metabolizers."
  (Source: `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/
  PA166104945>`__)

- Thioguanine

  "Consider an alternate agent or extreme dose reduction of thioguanine for
  patients who are TPMT or NUDT15 poor metabolizers. Start at 50-80% of
  target dose for patients who are TPMT or NUDT15 intermediate metabolizers."
  (Source: `PharmGKB <https://www.pharmgkb.org/guidelineAnnotation/
  PA166104965>`__)

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
     - Priority
   * - Normal Metabolizer
     - \*1/\*1
     - Normal/Routine/Low Risk
   * - Intermediate Metabolizer
     - \*1/\*6
     - Normal/Routine/Low Risk
   * - Poor Metabolizer
     - \*6/\*27
     - Abnormal/Priority/High Risk
   * - Indeterminate
     - \*28/\*80
     - None

Recommendations for UGT1A1
--------------------------

- Atazanavir

  "The CPIC dosing guideline recommends considering advising individuals who
  carry two decreased function UGT1A1 alleles about a substantial likelihood
  of developing jaundice, which may cause non-adherence. The dosing guideline
  recommends that alternative agents be considered if the risk of
  non-adherence due to jaundice is high. The risk of discontinuation is low
  and very low for individuals carrying one, or no decreased function UGT1A1
  alleles, respectively." (Source: `PharmGKB <https://www.pharmgkb.org/
  guidelineAnnotation/PA166128738>`__)

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

  * - SV Alleles
    - SV Name
    - Genotype
    - Reference
    - Gene Model
    - GRCh37
    - GRCh38
    - Data Type
    - Source
    - Coriell ID
    - Version
    - Description
  * -
    - Normal
    - \*1/\*2
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT1A4-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT1A4-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT1A4-3.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA11993
    - 0.9.0
    -
  * - \*S1
    - NoncodingDel1
    - \*1/\*S1
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT1A4-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT1A4-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT1A4-1.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA19908
    - 0.9.0
    -
  * - \*S1
    - NoncodingDel1Hom
    - \*S1/\*S1
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT1A4-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT1A4-5.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT1A4-5.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - HG03479
    - 0.16.0
    -
  * - \*S2
    - NoncodingDel2
    - \*1/\*S2
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT1A4-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT1A4-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT1A4-2.png>`
    - WGS
    -
    -
    - 0.10.0
    -
  * - \*S3
    - NoncodingDup1
    - \*1/\*S3
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT1A4-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT1A4-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT1A4-4.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - NA18632
    - 0.13.0
    -

PyPGx was recently applied to the entire high-coverage WGS dataset from 1KGP
(N=2,504). Click `here <https://github.com/sbslee/1kgp-pgx-paper/blob/main/
sv-tables/UGT1A4.md>`__ to see individual SV calls for UGT1A4, and
corresponding copy number profiles and allele fraction profiles.

UGT2B15
=======

SV summary for UGT2B15
----------------------

Below is comprehensive summary of SV described from real NGS studies:

.. list-table::
  :header-rows: 1

  * - SV Alleles
    - SV Name
    - Genotype
    - Reference
    - Gene Model
    - GRCh37
    - GRCh38
    - Data Type
    - Source
    - Coriell ID
    - Version
    - Description
  * -
    - Normal
    - \*1/\*2
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT2B15-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT2B15-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT2B15-2.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - HG00589
    - 0.4.1
    -
  * - \*S4
    - WholeDel1
    - \*2/\*S4
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT2B15-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT2B15-5.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT2B15-5.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - NA19024
    - 0.13.0
    -
  * -
    - WholeDel2
    - Indeterminate
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT2B15-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT2B15-6.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT2B15-6.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - NA19786
    - 0.14.0
    -
  * -
    - WholeDup1
    - Indeterminate
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT2B15-5.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT2B15-7.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT2B15-7.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - NA19776
    - 0.14.0
    -
  * - \*S1
    - PartialDel1
    - \*4/\*S1
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT2B15-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT2B15-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT2B15-1.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA11993
    - 0.4.1
    -
  * - \*S2
    - PartialDel2
    - \*2/\*S2
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT2B15-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT2B15-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT2B15-3.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - NA19160
    - 0.12.0
    -
  * - \*S3
    - PartialDel3
    - \*1/\*S3
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT2B15-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT2B15-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT2B15-4.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - NA19189
    - 0.13.0
    -
  * -
    - PartialDup1
    - Indeterminate
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT2B15-5.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT2B15-8.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT2B15-8.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - NA20821
    - 0.14.0
    -
  * -
    - PartialDup2
    - Indeterminate
    -
    -
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT2B15-9.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT2B15-9.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - HG03082
    - 0.16.0
    -

PyPGx was recently applied to the entire high-coverage WGS dataset from 1KGP
(N=2,504). Click `here <https://github.com/sbslee/1kgp-pgx-paper/blob/main/
sv-tables/UGT2B15.md>`__ to see individual SV calls for UGT2B15, and
corresponding copy number profiles and allele fraction profiles.

UGT2B17
=======

SV summary for UGT2B17
----------------------

This gene is known to have an extremely high rate of gene deletion
polymorphism in the population and thus requires SV analysis.

Below is comprehensive summary of SV described from real NGS studies:

.. list-table::
  :header-rows: 1

  * - SV Alleles
    - SV Name
    - Genotype
    - Reference
    - Gene Model
    - GRCh37
    - GRCh38
    - Data Type
    - Source
    - Coriell ID
    - Version
    - Description
  * -
    - Normal
    - \*1/\*1
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT2B17-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT2B17-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT2B17-3.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA19178
    - 0.4.1
    -
  * - \*2
    - WholeDel1
    - \*1/\*2
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT2B17-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT2B17-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT2B17-1.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA18855
    - 0.4.1
    -
  * - \*2
    - WholeDel1Hom
    - \*2/\*2
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT2B17-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT2B17-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT2B17-2.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA18617
    - 0.4.1
    -
  * - \*S2
    - PartialDel2
    - \*1/\*S2
    -
    -
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT2B17-8.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT2B17-8.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - HG03127
    - 0.16.0
    -
  * - \*S3
    - PartialDel3
    - \*1/\*S3
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT2B17-5.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT2B17-6.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT2B17-6.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - NA20886
    - 0.14.0
    -
  * - \*2, \*S1
    - WholeDel1+PartialDel1
    - \*2/\*S1
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT2B17-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT2B17-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT2B17-4.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - NA19160
    - 0.12.0
    -
  * - \*2, \*S2
    - WholeDel1+PartialDel2
    - \*2/\*S2
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT2B17-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT2B17-5.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT2B17-5.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - NA19189
    - 0.13.0
    -
  * - \*2, \*S3
    - WholeDel1+PartialDel3
    - \*2/\*S3
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT2B17-6.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT2B17-7.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT2B17-7.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - NA21090
    - 0.15.0
    -

PyPGx was recently applied to the entire high-coverage WGS dataset from 1KGP
(N=2,504). Click `here <https://github.com/sbslee/1kgp-pgx-paper/blob/main/
sv-tables/UGT2B17.md>`__ to see individual SV calls for UGT2B17, and
corresponding copy number profiles and allele fraction profiles.
