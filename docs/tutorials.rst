Tutorials
*********

GeT-RM WGS tutorial
===================

In this tutorial I'll walk you through PyPGx's genotype analysis using whole genome sequencing (WGS) data. By the end of this tutorial, you will have learned how to perform genotype analysis for genes with or without structural variation (SV), accordingly. I will also show how PyPGx can handle genomic data from two different Genome Reference Consortium Human (GRCh) builds: GRCh37 and GRCh38.

Before beginning this tutorial, create a new directory and change to that directory:

.. code-block:: text

  $ mkdir getrm-wgs-tutorial
  $ cd getrm-wgs-tutorial

The Centers for Disease Control and Prevention–based Genetic Testing Reference Materials Coordination Program (GeT-RM) has established `genomic DNA reference materials <https://www.cdc.gov/labquality/get-rm/inherited-genetic-diseases-pharmacogenetics/pharmacogenetics.html>`__  to help the genetic testing community obtain characterized reference materials. In particular, GeT-RM has made WGS data for 70 of reference samples publicly available for download and use from the `European Nucleotide Archive <https://www.ebi.ac.uk/ena/browser/view/PRJEB19931>`__. We will be using this WGS dataset throughout the tutorial.

Because downloading the entire WGS dataset is not feasible for most users due to its file size (i.e. a 30x WGS sample ≈ 90 GB), I have prepared input files ranging from 2 KB to 17.6 MB, for both GRCh37 and GRCh38. You can download those from:

.. code-block:: text

  $ wget https://raw.githubusercontent.com/sbslee/pypgx-data/main/getrm-wgs-tutorial/grch37-variants.vcf.gz
  $ wget https://raw.githubusercontent.com/sbslee/pypgx-data/main/getrm-wgs-tutorial/grch37-depth-of-coverage.zip
  $ wget https://raw.githubusercontent.com/sbslee/pypgx-data/main/getrm-wgs-tutorial/grch37-control-statistics-VDR.zip
  $ wget https://raw.githubusercontent.com/sbslee/pypgx-data/main/getrm-wgs-tutorial/grch38-variants.vcf.gz
  $ wget https://raw.githubusercontent.com/sbslee/pypgx-data/main/getrm-wgs-tutorial/grch38-depth-of-coverage.zip
  $ wget https://raw.githubusercontent.com/sbslee/pypgx-data/main/getrm-wgs-tutorial/grch38-control-statistics-VDR.zip

Genotyping genes with SV
------------------------

The first gene we are going to genotype is *CYP2D6* which has almost 150 star alleles including those with SV (e.g. gene deletions, duplications, and hybrids). To this end, we will run PyPGx's next-generation sequencing (NGS) pipeline:

.. code-block:: text

    $ pypgx run-ngs-pipeline \
    CYP2D6 \
    grch37-CYP2D6-pipeline \
    --variants grch37-variants.vcf.gz \
    --depth-of-coverage grch37-depth-of-coverage.zip \
    --control-statistics grch37-control-statistics-VDR.zip

Above will create a number of archive files:

.. code-block:: text

    Saved VcfFrame[Imported] to: grch37-CYP2D6-pipeline/imported-variants.zip
    Saved VcfFrame[Phased] to: grch37-CYP2D6-pipeline/phased-variants.zip
    Saved VcfFrame[Consolidated] to: grch37-CYP2D6-pipeline/consolidated-variants.zip
    Saved SampleTable[Alleles] to: grch37-CYP2D6-pipeline/alleles.zip
    Saved CovFrame[ReadDepth] to: grch37-CYP2D6-pipeline/read-depth.zip
    Saved CovFrame[CopyNumber] to: grch37-CYP2D6-pipeline/copy-number.zip
    Saved SampleTable[CNVCalls] to: grch37-CYP2D6-pipeline/cnv-calls.zip
    Saved SampleTable[Genotypes] to: grch37-CYP2D6-pipeline/genotypes.zip
    Saved SampleTable[Phenotypes] to: grch37-CYP2D6-pipeline/phenotypes.zip
    Saved SampleTable[Results] to: grch37-CYP2D6-pipeline/results.zip

Now that's it! You have successfully genotyped *CYP2D6* with WGS data.

Genotyping genes without SV
---------------------------

Next, let's run the same pipeline to genotype other gene called *CYP3A5*:

.. code-block:: text

    $ pypgx run-ngs-pipeline \
    CYP3A5 \
    grch37-CYP3A5-pipeline \
    --variants grch37-variants.vcf.gz \
    --depth-of-coverage grch37-depth-of-coverage.zip \
    --control-statistics grch37-control-statistics-VDR.zip

As before, above will create a number of archive files:

.. code-block:: text

    Saved VcfFrame[Imported] to: grch37-CYP3A5-pipeline/imported-variants.zip
    Saved VcfFrame[Phased] to: grch37-CYP3A5-pipeline/phased-variants.zip
    Saved VcfFrame[Consolidated] to: grch37-CYP3A5-pipeline/consolidated-variants.zip
    Saved SampleTable[Alleles] to: grch37-CYP3A5-pipeline/alleles.zip
    Saved SampleTable[Genotypes] to: grch37-CYP3A5-pipeline/genotypes.zip
    Saved SampleTable[Phenotypes] to: grch37-CYP3A5-pipeline/phenotypes.zip
    Saved SampleTable[Results] to: grch37-CYP3A5-pipeline/results.zip

For this run, you will notice two things:

1. There are less output files than there was with *CYP2D6*.
2. You will have received a warning from PyPGx that says: "UserWarning: The user provided CovFrame[DepthOfCoverage] even though the target gene does not have any star alleles defined by SVs. PyPGx will ignore it."

GRCh37 vs. GRCh38
-----------------

Thus far, we have only considered GRCh37 data. But we can also run the pipeline for GRCh38 data:

.. code-block:: text

    $ pypgx run-ngs-pipeline \
    CYP3A5 \
    grch38-CYP3A5-pipeline \
    --variants grch38-variants.vcf.gz \
    --assembly GRCh38

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
