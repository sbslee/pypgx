Tutorials
*********

GeT-RM tutorial
===============

In this tutorial I will show you how to genotype the *CYP2D6* gene with whole genome sequencing (WGS) data.

The Centers for Disease Control and Prevention (CDC)–based Genetic Testing Reference Materials Coordination Program (GeT-RM) has established `genomic DNA reference materials <https://www.cdc.gov/labquality/get-rm/inherited-genetic-diseases-pharmacogenetics/pharmacogenetics.html>`__  to help the genetic testing community obtain characterized reference materials. In particular, GeT-RM has made WGS data for 70 of reference samples publicly available for download and use from the `European Nucleotide Archive <https://www.ebi.ac.uk/ena/browser/view/PRJEB19931>`__. This is the WGS data we will be using for the tutorial.

Because downloading the entire WGS dataset is not feasible for most users due to its file size, I prepared input files:

.. code-block:: text

  $ mkdir getrm-tutorial
  $ cd getrm-tutorial
  $ wget https://raw.githubusercontent.com/sbslee/pypgx-data/main/getrm-tutorial/variants.vcf.gz
  $ wget https://raw.githubusercontent.com/sbslee/pypgx-data/main/getrm-tutorial/depth-of-coverage.tsv.gz
  $ wget https://raw.githubusercontent.com/sbslee/pypgx-data/main/getrm-tutorial/control-statistics-VDR.zip

Next, run the next-generation sequencing (NGS) pipeline:

.. code-block:: text

  $ pypgx run-ngs-pipeline \
  CYP2D6 \
  CYP2D6-pipeline \
  --variants variants.vcf.gz \
  --depth-of-coverage depth-of-coverage.tsv.gz \
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
