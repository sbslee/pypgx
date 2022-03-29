FAQ
***

This page provides answers to frequently asked questions (FAQ) for PyPGx.
Please open a GitHub issue `here <https://github.com/sbslee/pypgx/issues>`__
if you don't find answer to your question.

Reference vs. default star alleles
==================================

There are some confusions over reference vs. default star alleles, like this
GitHub issue (:issue:`51`). In the PGx field, genotypes and enzyme functions
are always reported relative to the gene's reference allele, instead of human
assembly such as GRCh37 or GRCh38. For example, the CYP2D6 gene has \*1 as
reference allele. Now, if you look at the CYP2D6 sequence of GRCh37, you will
find that it actually matches that of \*2. When you do the same for GRCh38,
its sequence matches that of \*1. Therefore, the default CYP2D6 alleles for
GRCh37 and GRCh38 are \*2 and \*1, respectively.

Flipped amino acid changes
==========================

Some users have asked why the PyPGx variant database sometimes has amino acid
changes that are opposite to what they see in other databases such as dbSNP.
For example, a GitHub issue (:issue:`52`) asks why PyPGx has the impact of
rs2032582 in the ABCB1 gene as Ala893Ser instead of Ser893Ala in dbSNP. This
is because the reference star allele, ABCB1\*1, has 893Ala while the ABCB1
sequence for both GRCh37 and GRCh38 has 893Ser (i.e. \*2). See
:ref:`faq:Reference vs. default star alleles` for more information.
