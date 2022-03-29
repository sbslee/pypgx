FAQ
***

This page provides answers to frequently asked questions (FAQ) for PyPGx.
Please open a GitHub issue `here <https://github.com/sbslee/pypgx/issues>`__
if you don't find answer to your question.

Reference vs. default star alleles
==================================

In the pharmacogenetic field, haplotypes and enzyme functions are always
reported relative to the gene's reference allele, instead of human assembly
such as GRCh37 or GRCh38. For example, the CYP2D6 gene has \*1 as reference
allele. Now, if you look at the CYP2D6 sequence of GRCh37, you will find that
it actually matches that of \*2. When you do the same for GRCh38, its
sequence matches that of \*1. Therefore, the default CYP2D6 alleles for
GRCh37 and GRCh38 are \*2 and \*1, respectively. This also means, in order
for a sample to have a \*1/\*1 diplotype in GRCh37, it will need to have some
homozygous variants first. This may sound a little bit strange at first (i.e.
reference star allele requiring variants), but I promise the more you think
about it, the better it will make sense. Related GitHub issues: :issue:`51`.

Flipped amino acid changes
==========================

In PyPGx you will notice that sometimes the variant database has amino acid
changes that are opposite to what you see in other databases such as dbSNP.
For example, a user has asked why PyPGx has the impact of rs2032582 from the
ABCB1 gene as Ala893Ser instead of Ser893Ala in dbSNP. This is because the
reference star allele, ABCB1\*1, has 893Ala while the ABCB1 sequence for both
GRCh37 and GRCh38 has 893Ser (i.e. \*2). See :ref:`faq:Reference vs. default
star alleles` for more information. Related GitHub issues: :issue:`52`.

Multiple alleles in single haplotype
====================================

It's possible, and fairly common, for single haplotype to fit the patterns of
multiple star alleles. For example, let's say there are two imaginary star
alleles \*A (defined by SNP1 and SNP2) and \*B (defined by SNP3). When the
three variants are haplotype phased together in cis, then the resulting
haplotype will have two candidate star alleles to choose from: \*A and \*B.
Now, PyPGx makes an opinionated pick to determine the final allele for
reporting purposes. Basically, PyPGx picks the allele with the highest
"priority" using the :meth:`pypgx.sort_alleles` method (`documentation
<https://pypgx.readthedocs.io/en/latest/api.html#pypgx.api.core.
sort_alleles>`__).
