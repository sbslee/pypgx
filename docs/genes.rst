Genes
*****

Introduction
============

This page describes gene-specific information. PyPGx currently supports
genotyping of a total of 58 pharmacogenes.

Many of the genes are known to have :ref:`structural variation (SV)
<glossary:Structural variation (SV)>` including gene deletions, duplications,
and hybrids. Please read the :ref:`readme:Structural variation detection`
page for more details.

Some genes have a genotype-phenotype table available from the
:ref:`glossary:Clinical Pharmacogenetics Implementation Consortium (CPIC)` or
the :ref:`glossary:Pharmacogenomics Knowledge Base (PharmGKB)`. Please read
the :ref:`readme:Phenotype prediction` page for more details.

Below is a summary table:

.. list-table::
   :header-rows: 1

   * - Gene
     - SV
     - Phenotype
     - PharmVar
     - CPIC
     - GRCh37
     - GRCh38
     - Notes
   * - ABCB1
     -
     -
     -
     -
     - `chr7:87130178-87345639 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A87130178%2D87345639&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr7:87500862-87716323 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A87500862%2D87716323&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CACNA1S`
     -
     - ✅
     -
     - ✅
     - `chr1:201005639-201084694 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A201005639%2D201084694&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr1:201036511-201115426 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A201036511%2D201115426&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CFTR`
     -
     - ✅
     -
     - ✅
     - `chr7:117117016-117311719 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A117117016%2D117311719&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr7:117477024-117671665 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A117477024%2D117671665&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP1A1`
     -
     -
     - ✅
     -
     - `chr15:75008882-75020951 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr15%3A75008882%2D75020951&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr15:74716541-74728528 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr15%3A74716541%2D74728528&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP1A2`
     -
     -
     - ✅
     -
     - `chr15:75038183-75051941 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr15%3A75038183%2D75051941&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr15:74745844-74759607 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr15%3A74745844%2D74759607&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP1B1`
     -
     -
     - ✅
     -
     - `chr2:38291745-38306323 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2%3A38291745%2D38306323&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr2:38064602-38079181 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2%3A38064602%2D38079181&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP2A6`
     - ✅
     -
     - ✅
     -
     - `chr19:41339442-41396352 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A41339442%2D41396352&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr19:40833540-40890447 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A40833540%2D40890447&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - CYP2A6 has pseudogene (CYP2A7).
   * - :ref:`genes:CYP2A13`
     -
     -
     - ✅
     -
     - `chr19:41574355-41622100 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A41574355%2D41622100&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr19:41068450-41116195 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A41068450%2D41116195&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP2B6`
     - ✅
     - ✅
     - ✅
     - ✅
     - `chr19:41427203-41534301 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A41427203%2D41534301&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr19:40921281-41028398 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A40921281%2D41028398&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - CYP2B6 has pseudogene (CYP2B7).
   * - :ref:`genes:CYP2C8`
     -
     -
     - ✅
     -
     - `chr10:96793528-96832254 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A96793528%2D96832254&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr10:95033771-95072497 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A95033771%2D95072497&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP2C9`
     -
     - ✅
     - ✅
     - ✅
     - `chr10:96695414-96752148 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A96695414%2D96752148&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr10:94935657-94993091 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A94935657%2D94993091&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP2C19`
     -
     - ✅
     - ✅
     - ✅
     - `chr10:96519437-96615962 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A96519437%2D96615962&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr10:94759680-94858547 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A94759680%2D94858547&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP2D6`
     - ✅
     - ✅
     - ✅
     - ✅
     - `chr22:42512500-42551883 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr22%3A42512500%2D42551883&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr22:42116498-42155810 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr22%3A42116498%2D42155810&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - CYP2D6 has pseudogene (CYP2D7).
   * - :ref:`genes:CYP2E1`
     - ✅
     -
     - ✅
     -
     - `chr10:135330866-135362620 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A135330866%2D135362620&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr10:133517362-133549123 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A133517362%2D133549123&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP2F1`
     -
     -
     - ✅
     -
     - `chr19:41617336-41637286 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A41617336%2D41637286&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr19:41111431-41131381 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A41111431%2D41131381&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP2J2`
     -
     -
     - ✅
     -
     - `chr1:60355979-60395470 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A60355979%2D60395470&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr1:59890307-59929773 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A59890307%2D59929773&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP2R1`
     -
     -
     - ✅
     -
     - `chr11:14896554-14916751 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr11%3A14896554%2D14916751&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr11:14875008-14895205 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr11%3A14875008%2D14895205&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP2S1`
     -
     -
     - ✅
     -
     - `chr19:41696111-41716444 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A41696111%2D41716444&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr19:41190218-41210539 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A41190218%2D41210539&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP2W1`
     -
     -
     - ✅
     -
     - `chr7:1019834-1032276 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A1019834%2D1032276&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr7:980180-992640 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A980180%2D992640&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP3A4`
     -
     -
     - ✅
     -
     - `chr7:99351582-99384811 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A99351582%2D99384811&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr7:99753966-99787184 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A99753966%2D99787184&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP3A5`
     -
     - ✅
     - ✅
     - ✅
     - `chr7:99242811-99280649 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A99242811%2D99280649&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr7:99645193-99682996 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A99645193%2D99682996&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP3A7`
     -
     -
     - ✅
     -
     - `chr7:99299659-99335823 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A99299659%2D99335823&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr7:99702035-99738196 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A99702035%2D99738196&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP3A43`
     -
     -
     - ✅
     -
     - `chr7:99422635-99466727 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A99422635%2D99466727&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr7:99825012-99869093 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A99825012%2D99869093&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP4A11`
     -
     -
     - ✅
     -
     - `chr1:47391859-47410148 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A47391859%2D47410148&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr1:46926187-46944476 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A46926187%2D46944476&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP4A22`
     -
     -
     - ✅
     -
     - `chr1:47600112-47618399 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A47600112%2D47618399&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr1:47134440-47152727 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A47134440%2D47152727&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP4B1`
     -
     -
     - ✅
     -
     - `chr1:47261669-47288021 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A47261669%2D47288021&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr1:46796045-46822413 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A46796045%2D46822413&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP4F2`
     - ✅
     -
     - ✅
     -
     - `chr19:15973833-16023930 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A15973833%2D16023930&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr19:15863022-15913074 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A15863022%2D15913074&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP17A1`
     -
     -
     - ✅
     -
     - `chr10:104587287-104600170 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A104587287%2D104600170&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr10:102827530-102840413 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A102827530%2D102840413&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP19A1`
     -
     -
     - ✅
     -
     - `chr15:51497253-51633795 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr15%3A51497253%2D51633795&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr15:51205056-51341596 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr15%3A51205056%2D51341596&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:CYP26A1`
     -
     -
     - ✅
     -
     - `chr10:94830646-94840641 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A94830646%2D94840641&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr10:93070892-93080885 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A93070892%2D93080885&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:DPYD`
     -
     - ✅
     - ✅
     - ✅
     - `chr1:97540298-98389615 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A97540298%2D98389615&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr1:97074742-97924034 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A97074742%2D97924034&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:F5`
     -
     - ✅
     -
     -
     - `chr1:169478188-169558719 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A169478188%2D169558719&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr1:169508950-169589481 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A169508950%2D169589481&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:G6PD`
     - ✅
     -
     -
     -
     - `chrX:153756604-153778233 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chrX%3A153756604%2D153778233&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chrX:154528389-154550018 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chrX%3A154528389%2D154550018&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - G6PD is located on X chromosome.
   * - :ref:`genes:GSTM1`
     - ✅
     -
     -
     -
     - `chr1:110227417-110239367 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A110227417%2D110239367&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr1:109684816-109696745 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A109684816%2D109696745&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - GSTP1
     -
     -
     -
     -
     - `chr11:67348065-67357124 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr11%3A67348065%2D67357124&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr11:67580811-67589653 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr11%3A67580811%2D67589653&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:GSTT1`
     - ✅
     -
     -
     -
     - `chr22:24373132-24387311 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr22%3A24373132%2D24387311&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr22_KI270879v1_alt:267307-281486 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr22_KI270879v1_alt%3A267307%2D281486&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - GSTT1 is located on different contigs between GRCh37 and GRCh38.
   * - :ref:`genes:IFNL3`
     -
     - ✅
     -
     -
     - `chr19:39731245-39738646 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A39731245%2D39738646&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr19:39240552-39248006 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A39240552%2D39248006&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - NAT1
     -
     -
     -
     -
     - `chr8:18064617-18084198 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr8%3A18064617%2D18084198&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr8:18207108-18226689 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr8%3A18207108%2D18226689&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - NAT2
     -
     -
     -
     -
     - `chr8:18245791-18261728 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr8%3A18245791%2D18261728&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr8:18388281-18404218 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr8%3A18388281%2D18404218&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:NUDT15`
     -
     - ✅
     - ✅
     - ✅
     - `chr13:48608702-48624364 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr13%3A48608702%2D48624364&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr13:48034725-48050221 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr13%3A48034725%2D48050221&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:POR`
     -
     -
     - ✅
     -
     - `chr7:75541419-75619173 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A75541419%2D75619173&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr7:75912154-75989855 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A75912154%2D75989855&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:PTGIS`
     -
     -
     - ✅
     -
     - `chr20:48117410-48187674 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr20%3A48117410%2D48187674&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr20:49500873-49571137 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr20%3A49500873%2D49571137&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:RYR1`
     -
     - ✅
     - ✅
     -
     - `chr19:38921339-39081204 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A38921339%2D39081204&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr19:38430690-38590564 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A38430690%2D38590564&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - SLC15A2
     -
     -
     -
     -
     - `chr3:121610170-121666034 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr3%3A121610170%2D121666034&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr3:121891400-121947188 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr3%3A121891400%2D121947188&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:SLC22A2`
     - ✅
     -
     -
     -
     - `chr6:160627786-160689853 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr6%3A160627786%2D160689853&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr6:160206754-160268821 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr6%3A160206754%2D160268821&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:SLCO1B1`
     -
     - ✅
     - ✅
     - ✅
     - `chr12:21281127-21395730 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr12%3A21281127%2D21395730&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr12:21128193-21242796 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr12%3A21128193%2D21242796&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - SLCO1B3
     -
     -
     -
     -
     - `chr12:20960637-21072845 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr12%3A20960637%2D21072845&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr12:20807704-20919911 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr12%3A20807704%2D20919911&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - SLCO2B1
     -
     -
     -
     -
     - `chr11:74859151-74920594 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr11%3A74859151%2D74920594&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr11:75148106-75209549 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr11%3A75148106%2D75209549&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:SULT1A1`
     - ✅
     -
     -
     -
     - `chr16:28601907-28636365 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr16%3A28601907%2D28636365&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr16:28590586-28625044 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr16%3A28590586%2D28625044&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:TBXAS1`
     -
     -
     - ✅
     -
     - `chr7:139525951-139723125 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A139525951%2D139723125&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr7:139826263-140023321 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A139826263%2D140023321&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:TPMT`
     -
     - ✅
     -
     - ✅
     - `chr6:18125541-18158400 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr6%3A18125541%2D18158400&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr6:18125310-18158169 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr6%3A18125310%2D18158169&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:UGT1A1`
     -
     - ✅
     -
     - ✅
     - `chr2:234662918-234687945 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2%3A234662918%2D234687945&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr2:233754269-233779300 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2%3A233754269%2D233779300&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:UGT1A4`
     - ✅
     -
     -
     -
     - `chr2:234624437-234684945 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2%3A234624437%2D234684945&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr2:233715735-233776300 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2%3A233715735%2D233776300&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - UGT2B7
     -
     -
     -
     -
     - `chr4:69959191-69981705 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr4%3A69959191%2D69981705&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr4:69093473-69115987 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr4%3A69093473%2D69115987&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:UGT2B15`
     - ✅
     -
     -
     -
     - `chr4:69506314-69542494 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr4%3A69506314%2D69542494&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr4:68640596-68676652 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr4%3A68640596%2D68676652&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - :ref:`genes:UGT2B17`
     - ✅
     -
     -
     -
     - `chr4:69399901-69437245 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr4%3A69399901%2D69437245&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr4:68534183-68571527 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr4%3A68534183%2D68571527&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - VKORC1
     -
     -
     -
     - ✅
     - `chr16:31099162-31109320 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr16%3A31099162%2D31109320&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr16:31087853-31097797 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr16%3A31087853%2D31097797&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     -
   * - XPC
     -
     -
     -
     -
     - `chr3:14183646-14223172 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr3%3A14183646%2D14223172&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
     - `chr3:14142146-14181672 <https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr3%3A14142146%2D14181672&hgsid=1251392659_FCwuNEZja7PPePnsIvfT1wF8Ke9Y>`__
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
     - Gene Model
     - GRCh37
     - GRCh38
     - Data Type
     - Source
     - Coriell ID
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
     -
   * - \*4
     - Deletion1Het
     - \*1/\*4
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-2.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-1.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-1.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA18617
     -
   * - \*4
     - Deletion1Hom
     - \*4/\*4
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-3.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-2.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-2.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA18952
     -
   * - \*4
     - Deletion2Het
     - \*1/\*4
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-2.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-6.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-6.png>`
     - WGS
     -
     -
     -
   * - \*4
     - Deletion3Het
     - \*4/\*9
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-2.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-7.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-7.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA18488
     -
   * - \*1x2
     - Duplication1
     - \*1x2/\*25
     - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-4.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-3.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-3.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA18861
     -
   * - \*1x2
     - Duplication2
     - \*1x2/\*2
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-4.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-10.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-10.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA12342
     -
   * - \*1x2
     - Duplication3
     - \*1x2/\*17
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-4.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-11.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-11.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA19129
     -
   * - Indeterminate
     - Hybrid1
     - Indeterminate
     -
     -
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-4.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-4.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - HG00436
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
     - \*12 has exons 1-2 of CYP2A7 origin and exons 3-9 of CYP2A6 origin (breakpoint in intron 2).
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
     - \*34 has axons 1-4 of CYP2A7 origin and axons 5-9 of CYP2A6 origin (breakpoint in intron 4)
   * -
     - PseudogeneDuplication
     - \*1/\*18
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2A6-7.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2A6-12.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2A6-12.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA06985
     -

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
     - Gene Model
     - GRCh37
     - GRCh38
     - Data Type
     - Source
     - Coriell ID
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
     -
   * - \*29
     - Hybrid
     - \*6/\*29
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2B6-2.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2B6-1.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2B6-1.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA19178
     - \*29 has exons 1-4 of CYP2B7 origin and exons 5-9 of CYP2A6 origin (breakpoint in intron 4).
   * - \*22x2
     - Duplication
     - \*6/\*22x2
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2B6-3.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2B6-3.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2B6-3.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA19190
     -

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
     - Gene Model
     - GRCh37
     - GRCh38
     - Data Type
     - Source
     - Coriell ID
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
     -
   * - \*5
     - DeletionHet
     - \*5/\*29
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2D6-2.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-1.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-1.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA18861
     -
   * - \*5
     - DeletionHom
     - \*5/\*5
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2D6-3.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-6.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-6.png>`
     - WGS
     -
     -
     -
   * - \*4x2
     - Duplication
     - \*2/\*4x2
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2D6-4.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-2.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-2.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA19819
     -
   * - \*1x3
     - Multiplication
     - \*1x3/\*10
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2D6-5.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-12.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-12.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA19190
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
     - \*68 has exon 1 of CYP2D6 origin and exons 2-9 of CYP2D7 origin (breakpoint in intron 1).
   * - \*68+\*4, \*68+\*4
     - Tandem1B
     - \*68+\*4/\*68+\*4
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2D6-7.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-13.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-13.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA12282
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
     -
   * - \*2x2, \*68+\*4
     - Duplication,Tandem1A
     - \*2x2/\*68+\*4
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2D6-12.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-10.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-10.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA21781
     -
   * - \*5, \*68+\*4
     - DeletionHet,Tandem1A
     - \*5/\*68+\*4
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2D6-13.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2D6-11.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2D6-11.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - HG01190
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
     - Gene Model
     - GRCh37
     - GRCh38
     - Data Type
     - Source
     - Coriell ID
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
     -
   * - \*S1
     - PartialDuplicationHet
     - \*1/\*S1
     - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2E1-2.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2E1-1.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2E1-1.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA19920
     - \*S1 is linked to \*7.
   * - \*S1
     - PartialDuplicationHom
     - \*S1/\*S1
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2E1-5.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2E1-7.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2E1-7.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA19309
     -
   * - \*1x2
     - Duplication1
     - \*1/\*1x2
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2E1-3.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2E1-4.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2E1-4.png>`
     - WGS
     -
     -
     -
   * - \*7x2
     - Duplication1
     - \*1/\*7x2
     - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2E1-3.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2E1-2.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2E1-2.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA19095
     -
   * - \*1x2
     - Duplication2
     - \*1x2/\*7
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2E1-3.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2E1-6.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2E1-6.png>`
     - WGS
     - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
     - NA19225
     -
   * - \*7x3
     - Multiplication
     - \*7/\*7x3
     -
     - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP2E1-4.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP2E1-3.png>`
     - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP2E1-3.png>`
     - WGS
     - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
     - NA19908
     -

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
    - Gene Model
    - GRCh37
    - GRCh38
    - Data Type
    - Source
    - Coriell ID
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
    -
  * - \*DEL
    - DeletionHet
    - \*1/\*DEL
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-CYP4F2-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-CYP4F2-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-CYP4F2-1.png>`
    - WGS
    -
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

  * - Star Allele
    - SV Name
    - Genotype
    - Reference
    - Gene Model
    - GRCh37
    - GRCh38
    - Data Type
    - Source
    - Coriell ID
    - Description
  * -
    - Female
    - \*B/\*B
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-G6PD-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-G6PD-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-G6PD-1.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - HG00276
    -
  * - \*MALE
    - Male
    - \*B/\*MALE
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-G6PD-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-G6PD-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-G6PD-2.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - HG00436
    -

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
    - Gene Model
    - GRCh37
    - GRCh38
    - Data Type
    - Source
    - Coriell ID
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
    -
  * - \*0
    - DeletionHet
    - \*0/\*A
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-GSTM1-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTM1-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTM1-1.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA18855
    -
  * - \*0
    - DeletionHom
    - \*0/\*0
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-GSTM1-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTM1-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTM1-2.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA10831
    -
  * - \*Ax2
    - Duplication
    - \*A/\*Ax2
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-GSTM1-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTM1-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTM1-3.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA19908
    -
  * - \*Bx2
    - Duplication
    - \*A/\*Bx2
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-GSTM1-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTM1-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTM1-4.png>`
    - WGS
    -
    -
    -
  * -
    - UpstreamDeletionHet
    - \*A/\*B
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-GSTM1-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTM1-6.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTM1-6.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - NA19005
    -
  * - \*0
    - DeletionHet,UpstreamDeletionHet
    - \*0/\*A
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-GSTM1-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTM1-7.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTM1-7.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - NA06984
    -

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

  * - Star Allele
    - SV Name
    - Genotype
    - Reference
    - Gene Model
    - GRCh37
    - GRCh38
    - Data Type
    - Source
    - Coriell ID
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
    -
  * - \*0
    - DeletionHet
    - \*0/\*A
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-GSTT1-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTT1-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTT1-1.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA19908
    -
  * - \*0
    - DeletionHom
    - \*0/\*0
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-GSTT1-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-GSTT1-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-GSTT1-2.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA11832
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
    - Gene Model
    - GRCh37
    - GRCh38
    - Data Type
    - Source
    - Coriell ID
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
    -
  * - \*S1
    - Intron9Deletion
    - \*1/\*S1
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-SLC22A2-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SLC22A2-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SLC22A2-2.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA18855
    -
  * - \*S2
    - Exon11Deletion
    - \*1/\*S2
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-SLC22A2-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SLC22A2-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SLC22A2-3.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA19819
    -
  * - \*S1, \*S2
    - Intron9Deletion,Exon11Deletion
    - \*S1/\*S2
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-SLC22A2-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SLC22A2-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SLC22A2-4.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - NA19030
    -

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
    - Gene Model
    - GRCh37
    - GRCh38
    - Data Type
    - Source
    - Coriell ID
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
    -
  * - \*DEL
    - DeletionHet
    - \*1/\*DEL
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-SULT1A1-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SULT1A1-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SULT1A1-1.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA18942
    -
  * - \*1x2
    - Duplication
    - \*1x2/\*2
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-SULT1A1-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SULT1A1-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SULT1A1-2.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA18509
    -
  * - \*1x3
    - Multiplication1
    - \*1x3/\*2
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-SULT1A1-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SULT1A1-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SULT1A1-3.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA18868
    -
  * - \*1x4
    - Multiplication2
    - \*1x4/\*2
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-SULT1A1-5.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SULT1A1-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SULT1A1-4.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA18484
    -
  * - \*1x3, \*2x2
    - Multiplication2
    - \*1x3/\*2x2
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-SULT1A1-6.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-SULT1A1-5.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-SULT1A1-5.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA19143
    -

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
    - Gene Model
    - GRCh37
    - GRCh38
    - Data Type
    - Source
    - Coriell ID
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
    -
  * - \*S1
    - Intron1DeletionA
    - \*1/\*S1
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT1A4-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT1A4-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT1A4-1.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA19908
    -
  * - \*S2
    - Intron1DeletionB
    - \*1/\*S2
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT1A4-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT1A4-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT1A4-2.png>`
    - WGS
    -
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
    - Gene Model
    - GRCh37
    - GRCh38
    - Data Type
    - Source
    - Coriell ID
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
    -
  * - \*S1
    - PartialDeletion1
    - \*4/\*S1
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT2B15-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT2B15-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT2B15-1.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA11993
    -
  * - \*S2
    - PartialDeletion2
    - \*2/\*S2
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT2B15-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT2B15-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT2B15-3.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - NA19160
    -

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
    - Gene Model
    - GRCh37
    - GRCh38
    - Data Type
    - Source
    - Coriell ID
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
    -
  * - \*2
    - DeletionHet
    - \*1/\*2
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT2B17-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT2B17-1.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT2B17-1.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA18855
    -
  * - \*2
    - DeletionHom
    - \*2/\*2
    - `Lee et al., 2019 <https://pubmed.ncbi.nlm.nih.gov/31206625/>`__
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT2B17-3.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT2B17-2.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT2B17-2.png>`
    - WGS
    - `GeT-RM <https://pubmed.ncbi.nlm.nih.gov/26621101/>`__
    - NA18617
    -
  * - \*2, \*S1
    - PartialDeletionHet
    - \*2/\*S1
    -
    - :download:`Model <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/gene-model-UGT2B17-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh37-UGT2B17-4.png>`
    - :download:`Profile <https://raw.githubusercontent.com/sbslee/pypgx-data/main/dpsv/GRCh38-UGT2B17-4.png>`
    - WGS
    - `1KGP <https://www.biorxiv.org/content/10.1101/2021.02.06.430068v2>`__
    - NA19160
    -
