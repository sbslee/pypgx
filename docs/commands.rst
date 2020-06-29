Commands
********

pgkb
====

Synopsis
--------

pypgx pgkb [options]

Description
-----------

Extract CPIC guidelines using PharmGKB API.

Options
-------

-h, --help  show this help message and exit
-o FILE     output to FILE [stdout]
-t          will only look first three guidelines

report
======

Synopsis
--------

pypgx report [options] gt

Description
-----------

Create HTML report using data from Stargazer. ``gt`` is input genotype file 
from Stargazer.

Options
-------

-h, --help  show this help message and exit
-o FILE     output to FILE [stdout]

sdf2gdf
=======

Synopsis
--------

pypgx sdf2gdf [options] sdf id [id ...]

Description
-----------

Create GDF file from SDF file. ``sdf`` is input SDF file. ``id`` is sample ID.

Options
-------

-h, --help  show this help message and exit
-o FILE     output to FILE [stdout]

bam2sdf
=======

Synopsis
--------

pypgx bam2sdf [options] tg cg bam [bam ...]

Description
-----------

Create SDF file from BAM file(s). ``tg`` is target gene. ``cg`` is control 
gene. ``bam`` is input BAM file.

Options
-------

-h, --help  show this help message and exit
-o FILE     output to FILE [stdout]

bam2gdf
=======

Synopsis
--------

pypgx bam2gdf [options] tg cg bam [bam ...]

Description
-----------

Create GDF file from BAM file(s). ``tg`` is target gene. ``cg`` is control 
gene. ``bam`` is input BAM file.

Options
-------

-h, --help  show this help message and exit
-o FILE     output to FILE [stdout]

minivcf
=======

Synopsis
--------

pypgx minivcf [options] vcf region

Description
-----------

Slice VCF file. ``vcf`` is input VCF file. ``region`` is target region.

Options
-------

-h, --help  show this help message and exit
-o FILE     output to FILE [stdout]

merge
========

Synopsis
--------

pypgx merge [options] vcf [vcf ...]

Description
-----------

Merge VCF files. ``vcf`` is input VCF file.

Options
-------

-h, --help  show this help message and exit
-r STR      genomic region (chr:start-end)
-o FILE     output to FILE [stdout]