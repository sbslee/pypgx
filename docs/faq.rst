FAQ
***

This page provides answers to frequently asked questions (FAQ) for PyPGx.
Please open a GitHub issue `here <https://github.com/sbslee/pypgx/issues>`__
if you don't find answer to your question.

Reference vs. default star alleles
==================================

In the pharmacogenetic field, haplotypes and enzyme functions are always
reported relative to the gene's reference star allele (some people refer it
as "wild-type" allele, but reference allele is the preferred term), instead
of human assembly such as GRCh37 or GRCh38. For example, the CYP2D6 gene has
CYP2D6\*1 as reference allele. Now, if you look at the CYP2D6 sequence of
GRCh37, you will find that it actually matches that of CYP2D6\*2. When you do
the same for GRCh38, its sequence matches that of CYP2D6\*1. Therefore, the
"default" alleles for GRCh37 and GRCh38 are CYP2D6\*2 and \*1, respectively.
This also means, in order for a sample to have a CYP2D6\*1/\*1 diplotype in
GRCh37, it will need to present some homozygous variants. This may sound a
little bit strange at first (i.e. reference star allele requiring variants to
be called), but I promise the more you think about it, the better it will
make sense. Related GitHub issues: :issue:`51`.

Flipped amino acid changes
==========================

In PyPGx you will notice that sometimes the variant database has amino acid
changes that are opposite to what you see in other databases such as dbSNP.
For example, PyPGx reports that the impact of rs2032582 from the ABCB1 gene
as Ala893Ser instead of Ser893Ala in dbSNP. This is because the reference
star allele, ABCB1\*1, has 893Ala while the ABCB1 sequence for both GRCh37
and GRCh38 has 893Ser (i.e. ABCB1\*2). See :ref:`faq:Reference vs. default
star alleles` for more information. Related GitHub issues: :issue:`52`.

Multiple alleles in single haplotype
====================================

It's possible, and fairly common, for single haplotype to fit the patterns of
multiple star alleles. For example, let's say there are two imaginary star
alleles of your favorite gene: \*A (which is defined by SNP1 and SNP2) and
\*B (which is defined by SNP3). When the three variants are haplotype phased
together in cis, then the resulting haplotype will have two candidate star
alleles to choose from: \*A and \*B. Now, PyPGx makes an opinionated pick to
determine the final allele for reporting purposes. Basically, PyPGx picks the
allele with the highest "priority" using the :meth:`pypgx.sort_alleles`
method (`documentation <https://pypgx.readthedocs.io/en/latest/api.html#pypgx
.api.core.sort_alleles>`__).

Variant caller choice
=====================

Starting with the 0.14.0 version, it is strongly recommended to use the
:command:`create-input-vcf` command for creating a VCF file (containing
SNVs/indels) from BAM files before running the `NGS pipeline <https://pypgx.
readthedocs.io/en/latest/readme.html#ngs-pipeline>`__. Prior to this command,
users had been instructed to create input VCF file from BAM files on their
own using a variant caller of their choice (e.g. GATK4, bcftools, DRAGEN,
DeepVariant). This can raise several potential problems such as decreased
reproducibility of PyPGx results and users providing incorrectly formatted
VCF to PyPGx. Another problem is the assumption that users are already
familiar with the ins and outs of variant calling (e.g. know how to control
the balance between sensitivity vs. specificity). PyPGx strongly recommends
producing input VCF that contains all possible SNVs/indels to achieve maximum
sensitivity because it will only use known variants for genotyping purposes
anyways (i.e. variants used to define star alleles); therefore, PyPGx is
actually quite robust against false positives. By introducing
:command:`create-input-vcf`, the main goal was to help standardize the NGS
pipeline even further. Related GitHub issues: :issue:`54`.

Disclaimer: There will be cases where more sophisticated variant callers are
preferred, or even suitable, to generate input VCF for PyPGx. For example, if
your samples are ancient DNA then you would probably want to use callers like
ATLAS to correct for ancient DNA damage. Also, if you've already created
input VCF for purposes other than running PyPGx and you want/need to be
consistent with the other variant-level analyses you may also just use the
same VCF for PyPGx. The bottom line is, if you are going to create your own
input VCF, then you need to know what you are doing. Otherwise, it's probably
safer to use :command:`create-input-vcf`.

``chr22_KI270879v1_alt`` in GRCh38
==================================

Users may encounter an error like below when working with GRCh38 data:

.. code-block:: text

    $ pypgx prepare-depth-of-coverage \
    depth-of-coverage.zip \
    HG00276_PyPGx.sorted.markdup.recal.bam \
    --assembly GRCh38
    Traceback (most recent call last):
      File "/Users/sbslee/opt/anaconda3/envs/fuc/bin/pypgx", line 33, in <module>
        sys.exit(load_entry_point('pypgx', 'console_scripts', 'pypgx')())
      File "/Users/sbslee/Desktop/pypgx/pypgx/__main__.py", line 33, in main
        commands[args.command].main(args)
      File "/Users/sbslee/Desktop/pypgx/pypgx/cli/prepare_depth_of_coverage.py", line 90, in main
        archive = utils.prepare_depth_of_coverage(
      File "/Users/sbslee/Desktop/pypgx/pypgx/api/utils.py", line 1247, in prepare_depth_of_coverage
        cf = pycov.CovFrame.from_bam(bams, regions=regions, zero=True)
      File "/Users/sbslee/Desktop/fuc/fuc/api/pycov.py", line 345, in from_bam
        results += pysam.depth(*(bams + args + ['-r', region]))
      File "/Users/sbslee/opt/anaconda3/envs/fuc/lib/python3.9/site-packages/pysam/utils.py", line 69, in __call__
        raise SamtoolsError(
    pysam.utils.SamtoolsError: 'samtools returned with error 1: stdout=, stderr=samtools depth: cannot parse region "chr22_KI270879v1_alt:267307-281486"\n'

This is a GRCh38-specific issue. One of the genes with SV is GSTT1 and it is
located in the contig ``chr22_KI270879v1_alt``, which is missing in input BAM
file. That's why the :command:`prepare-depth-of-coverage` command is
complaining. For more details, please see the following articles:
:ref:`readme:GRCh37 vs. GRCh38` and :ref:`genes:GRCh38 data for GSTT1`.
Related GitHub issues: :issue:`65`.
