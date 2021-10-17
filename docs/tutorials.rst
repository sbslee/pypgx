Tutorials
*********

GeT-RM WGS tutorial
===================

In this tutorial I will show you how to genotype the *CYP2D6* gene with whole genome sequencing (WGS) data.

The Centers for Disease Control and Prevention (CDC)–based Genetic Testing Reference Materials Coordination Program (GeT-RM) has established `genomic DNA reference materials <https://www.cdc.gov/labquality/get-rm/inherited-genetic-diseases-pharmacogenetics/pharmacogenetics.html>`__  to help the genetic testing community obtain characterized reference materials. In particular, GeT-RM has made WGS data for 70 of reference samples publicly available for download and use from the `European Nucleotide Archive <https://www.ebi.ac.uk/ena/browser/view/PRJEB19931>`__. This is the WGS data we will be using for the tutorial.

Because downloading the entire WGS dataset is not feasible for most users due to its file size, I prepared input files:

.. code-block:: text

  $ mkdir getrm-wgs-tutorial
  $ cd getrm-wgs-tutorial
  $ wget https://raw.githubusercontent.com/sbslee/pypgx-data/main/getrm-wgs-tutorial/variants.vcf.gz
  $ wget https://raw.githubusercontent.com/sbslee/pypgx-data/main/getrm-wgs-tutorial/depth-of-coverage.zip
  $ wget https://raw.githubusercontent.com/sbslee/pypgx-data/main/getrm-wgs-tutorial/control-statistics-VDR.zip

Next, run the next-generation sequencing (NGS) pipeline:

.. code-block:: text

  $ pypgx run-ngs-pipeline \
  CYP2D6 \
  CYP2D6-pipeline \
  --variants variants.vcf.gz \
  --depth-of-coverage depth-of-coverage.zip \
  --control-statistics control-statistics-VDR.zip

Above will create a number of archive files:

.. code-block:: text

  Saved VcfFrame[Imported] to: CYP2D6-pipeline/imported-variants.zip
  Saved VcfFrame[Phased] to: CYP2D6-pipeline/phased-variants.zip
  Saved VcfFrame[Consolidated] to: CYP2D6-pipeline/consolidated-variants.zip
  Saved SampleTable[Alleles] to: CYP2D6-pipeline/alleles.zip
  Saved CovFrame[ReadDepth] to: CYP2D6-pipeline/read-depth.zip
  Saved CovFrame[CopyNumber] to: CYP2D6-pipeline/copy-number.zip
  Saved SampleTable[CNVCalls] to: CYP2D6-pipeline/cnv-calls.zip
  Saved SampleTable[Genotypes] to: CYP2D6-pipeline/genotypes.zip
  Saved SampleTable[Results] to: CYP2D6-pipeline/results.zip

Now that's it! You have successfully genotyped *CYP2D6* with WGS data.

Coriell Affy tutorial
=====================

In this tutorial I will show you how to genotype the *CYP3A5* gene with chip data.

Coriell Institute has carried out Affy 6.0 genotyping on many of the 1000 Genomes Project (1KGP) samples whose data are available on 1KGP's `FTP site <http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/>`__. For this tutorial we will be using the file ``ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_no_ped.vcf.gz`` which contains variant data for 355 samples.

For convenience, I prepared input files:

.. code-block:: text

  $ mkdir coriell-affy-tutorial
  $ cd coriell-affy-tutorial
  $ wget https://raw.githubusercontent.com/sbslee/pypgx-data/main/coriell-affy-tutorial/variants.vcf.gz

Next, run the chip pipeline:

.. code-block:: text

  $ pypgx run-chip-pipeline \
  CYP3A5 \
  CYP3A5-pipeline \
  variants.vcf.gz

Above will create a number of archive files:

.. code-block:: text

  Saved VcfFrame[Imported] to: CYP3A5-pipeline/imported-variants.zip
  Saved VcfFrame[Phased] to: CYP3A5-pipeline/phased-variants.zip
  Saved VcfFrame[Consolidated] to: CYP3A5-pipeline/consolidated-variants.zip
  Saved SampleTable[Alleles] to: CYP3A5-pipeline/alleles.zip
  Saved SampleTable[Genotypes] to: CYP3A5-pipeline/genotypes.zip
  Saved SampleTable[Phenotypes] to: CYP3A5-pipeline/phenotypes.zip
  Saved SampleTable[Results] to: CYP3A5-pipeline/results.zip

Now that’s it! You have successfully genotyped *CYP3A5* with chip data.
