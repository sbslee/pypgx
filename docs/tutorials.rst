Tutorials
*********

This page provides useful tutorials for running PyPGx. Throughout the page,
it's assumed that you have installed the latest version of PyPGx and also
downloaded the appropriate resource bundle (i.e. matching version). For more
details, please see :ref:`readme:Installation` and :ref:`readme:Resource
bundle`.

GeT-RM WGS tutorial
===================

In this tutorial I'll walk you through PyPGx's genotype analysis using whole
genome sequencing (WGS) data. By the end of this tutorial, you will have
learned how to perform genotype analysis for genes with or without structural
variation (SV), accordingly. I will also show how PyPGx can handle genomic
data from two different Genome Reference Consortium Human (GRCh) builds:
GRCh37 (hg19) and GRCh38 (hg38).

Before beginning this tutorial, create a new directory and change to that
directory:

.. code-block:: text

    $ mkdir getrm-wgs-tutorial
    $ cd getrm-wgs-tutorial

The Centers for Disease Control and Prevention–based Genetic Testing
Reference Materials Coordination Program (GeT-RM) has established `genomic
DNA reference materials <https://www.cdc.gov/labquality/get-rm/inherited-
genetic-diseases-pharmacogenetics/pharmacogenetics.html>`__  to help the
genetic testing community obtain characterized reference materials. In
particular, GeT-RM has made WGS data for 70 of reference samples publicly
available for download and use from the `European Nucleotide Archive
<https://www.ebi.ac.uk/ena/browser/view/PRJEB19931>`__. We will be using this
WGS dataset throughout the tutorial.

Obtaining input files
---------------------

Because downloading the entire WGS dataset is probably not feasible for most
users due to large file size (i.e. a 30x WGS sample ≈ 90 GB), I have prepared
input files ranging from 2 KB to 25.5 MB, for both GRCh37 and GRCh38. You can
easily download these with:

.. code-block:: text

    $ wget https://raw.githubusercontent.com/sbslee/pypgx-data/main/getrm-wgs-tutorial/grch37-variants.vcf.gz
    $ wget https://raw.githubusercontent.com/sbslee/pypgx-data/main/getrm-wgs-tutorial/grch37-variants.vcf.gz.tbi
    $ wget https://raw.githubusercontent.com/sbslee/pypgx-data/main/getrm-wgs-tutorial/grch37-depth-of-coverage.zip
    $ wget https://raw.githubusercontent.com/sbslee/pypgx-data/main/getrm-wgs-tutorial/grch37-control-statistics-VDR.zip
    $ wget https://raw.githubusercontent.com/sbslee/pypgx-data/main/getrm-wgs-tutorial/grch38-variants.vcf.gz
    $ wget https://raw.githubusercontent.com/sbslee/pypgx-data/main/getrm-wgs-tutorial/grch38-variants.vcf.gz.tbi
    $ wget https://raw.githubusercontent.com/sbslee/pypgx-data/main/getrm-wgs-tutorial/grch38-depth-of-coverage.zip
    $ wget https://raw.githubusercontent.com/sbslee/pypgx-data/main/getrm-wgs-tutorial/grch38-control-statistics-VDR.zip

Let's look at the metadata for some of these files:

.. code-block:: text

    $ pypgx print-metadata grch37-depth-of-coverage.zip
    Assembly=GRCh37
    SemanticType=CovFrame[DepthOfCoverage]
    Platform=WGS

.. code-block:: text

    $ pypgx print-metadata grch38-control-statistics-VDR.zip
    Control=VDR
    Assembly=GRCh38
    SemanticType=SampleTable[Statistics]
    Platform=WGS

At this point, you are now ready to move on to the next step.

Optionally, in case you are interested in creating above input files on your
own, I have also prepared "mini" BAM files for GRCh37 where the original
sequencing data from GeT-RM have been sliced to contain genes used by PyPGx
only. You can download them `here <https://1drv.ms/u/
s!Apgoq3uQ2gCqgrovIFKJSi-ECXY9pw?e=uP5EeU>`__. You will also need reference
FASTA when creating input VCF, which can be downloaded from `here
<https://1drv.ms/u/s!Apgoq3uQ2gCqgt4qGq9YsumpVk9xJQ?e=ZewLHu>`__.

Once you are finished downloading the mini BAM files and the reference FASTA
file, let's first create input VCF:

.. code-block:: text

    $ pypgx create-input-vcf \
    grch37-variants.vcf.gz \
    /path/to/genome.fa \
    grch37-bam.list

Note that this step can take some time to run. For example, it takes about 1
hour to finish using my personal MacBook Air (M1, 2020) with 8 GB of memory.

Next, we will compute depth of coverage for genes that are known to have SV:

.. code-block:: text

    $ pypgx prepare-depth-of-coverage \
    grch37-depth-of-coverage.zip \
    grch37-bam.list

This step should be quick. It finishes in less than 30 seconds with my laptop.

Finally, we will compute control statistics using the VDR gene as control
locus, which is required when converting read depth to copy number:

.. code-block:: text

    $ pypgx compute-control-statistics \
    VDR \
    grch37-control-statistics-VDR.zip \
    grch37-bam.list

This step should be quick as well. It finishes in less than 5 seconds with my
laptop.

Genotyping genes with SV
------------------------

The first gene we are going to genotype is CYP2D6, which has almost 150
star alleles including those with SV (e.g. gene deletions, duplications, and
hybrids). To this end, we will run PyPGx's next-generation sequencing (NGS)
pipeline (see :ref:`readme:NGS pipeline` for more details):

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

In addition to these files, PyPGx will have also created two directories
called ``copy-number-profile`` and ``allele-fraction-profile``.

Now let's make sure the genotype results are correct by comparing them with the validation data:

.. code-block:: text

    $ wget https://raw.githubusercontent.com/sbslee/pypgx-data/main/getrm-wgs-tutorial/grch37-CYP2D6-results.zip
    $ pypgx compare-genotypes grch37-CYP2D6-pipeline/results.zip grch37-CYP2D6-results.zip
    # Genotype
    Total: 70
    Compared: 70
    Concordance: 1.000 (70/70)
    # CNV
    Total: 70
    Compared: 70
    Concordance: 1.000 (70/70)

That's it, you have successfully genotyped CYP2D6 with WGS data!

Genotyping genes without SV
---------------------------

The next gene we're going to genotype is CYP3A5. Unlike CYP2D6, this gene
does not have any star alleles with SV. Therefore, we only need to provide
``grch37-variants.vcf.gz`` to the NGS pipeline:

.. code-block:: text

    $ pypgx run-ngs-pipeline \
    CYP3A5 \
    grch37-CYP3A5-pipeline \
    --variants grch37-variants.vcf.gz

Above will create a number of archive files:

.. code-block:: text

    Saved VcfFrame[Imported] to: grch37-CYP3A5-pipeline/imported-variants.zip
    Saved VcfFrame[Phased] to: grch37-CYP3A5-pipeline/phased-variants.zip
    Saved VcfFrame[Consolidated] to: grch37-CYP3A5-pipeline/consolidated-variants.zip
    Saved SampleTable[Alleles] to: grch37-CYP3A5-pipeline/alleles.zip
    Saved SampleTable[Genotypes] to: grch37-CYP3A5-pipeline/genotypes.zip
    Saved SampleTable[Phenotypes] to: grch37-CYP3A5-pipeline/phenotypes.zip
    Saved SampleTable[Results] to: grch37-CYP3A5-pipeline/results.zip

Plus the ``allele-fraction-profile`` directory.

Now you have successfully genotyped CYP3A5 as well!

.. note::
    Note that if you provide ``grch37-depth-of-coverage.zip`` and
    ``grch37-control-statistics-VDR.zip`` to the pipeline, PyPGx will still
    run without any issues, but it will output a warning that says those
    files will be ignored. This is so that users don't have to memorize which
    gene requires SV analysis. In other words, users can provide the same
    input files for all target genes.

Genotyping with GRCh38 data
---------------------------

Thus far, we have only considered GRCh37 data. But we can also run the
pipeline for GRCh38 data by changing the ``--assembly`` option:

.. code-block:: text

    $ pypgx run-ngs-pipeline \
    CYP3A5 \
    grch38-CYP3A5-pipeline \
    --variants grch38-variants.vcf.gz \
    --assembly GRCh38

Which will create:

.. code-block:: text

    Saved VcfFrame[Imported] to: grch38-CYP3A5-pipeline/imported-variants.zip
    Saved VcfFrame[Phased] to: grch38-CYP3A5-pipeline/phased-variants.zip
    Saved VcfFrame[Consolidated] to: grch38-CYP3A5-pipeline/consolidated-variants.zip
    Saved SampleTable[Alleles] to: grch38-CYP3A5-pipeline/alleles.zip
    Saved SampleTable[Genotypes] to: grch38-CYP3A5-pipeline/genotypes.zip
    Saved SampleTable[Phenotypes] to: grch38-CYP3A5-pipeline/phenotypes.zip
    Saved SampleTable[Results] to: grch38-CYP3A5-pipeline/results.zip

Now let’s make sure the genotype results are correct by comparing them with
the GRCh37 results:

.. code-block:: text

    $ pypgx compare-genotypes grch37-CYP3A5-pipeline/results.zip grch38-CYP3A5-pipeline/results.zip
    # Genotype
    Total: 70
    Compared: 70
    Concordance: 1.000 (70/70)
    # CNV
    Total: 70
    Compared: 0
    Concordance: N/A

Congratulations, you have completed this tutorial!

Coriell Affy tutorial
=====================

In this tutorial I will show you how to genotype the CYP3A5 gene with chip data.

Coriell Institute has carried out Affy 6.0 genotyping on many of the 1000 Genomes Project (1KGP) samples whose data are available on 1KGP's `FTP site <http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/>`__. For this tutorial we will be using the file ``ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_no_ped.vcf.gz`` which contains variant data for 355 samples.

For convenience, I prepared input files:

.. code-block:: text

  $ mkdir coriell-affy-tutorial
  $ cd coriell-affy-tutorial
  $ wget https://raw.githubusercontent.com/sbslee/pypgx-data/main/coriell-affy-tutorial/variants.vcf.gz
  $ wget https://raw.githubusercontent.com/sbslee/pypgx-data/main/coriell-affy-tutorial/variants.vcf.gz.tbi

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

Now that’s it! You have successfully genotyped CYP3A5 with chip data.
