Changelog
*********

0.3.0
-----

* Updated the ``phenotyper()`` method for the CYP2C9 gene (`#18 <https://github.com/sbslee/pypgx/issues/18>`_).
* Changed the function of CYP2C9*3 in ``star_table.tsv`` from decreased function to no function.
* Updated the ``calculate-read-depth`` command to accept single BAM file as input (`#19 <https://github.com/sbslee/pypgx/issues/19>`_).

0.2.0
-----

* From now on, pypgx will be hosted on Anaconda instead of PyPI. This is because pypgx needs access to many non-Python programs and these are not available via PyPI.
* Updated the license from GPL-3.0 License to MIT License.
* Temporarily removed all the of commands except ``bam2vcf2``, ``bam2gdf``, and ``compgt``. Some of the removed commands will be brought back in the later versions.
* The ``bam2vcf2`` command has been renamed as ``call-variants-gatk-sge``.
* The ``bam2gdf`` command has been renamed as ``calculate-read-depth``.
* The ``compgt`` command has been renamed as ``compare-stargazer-calls``.

0.1.37
------

* Removed dependency to the ``vcfgo`` package.
* Removed the ``mergevcf`` command.

0.1.36
------

* Updated NUDT15 phenotype (`#9 <https://github.com/sbslee/pypgx/pull/9>`_).

0.1.35
------

* Updated phenotypes for TPMT, UGT1A1, and DPYD (`#5 <https://github.com/sbslee/pypgx/issues/5>`_).

0.1.34
------

* Added a new configuration parameter called `conda_env`.
* Added a conda recipe (conda/stargazer-conda.yml) for Stargazer environment.

0.1.33
------

* No significant changes.

0.1.32
------

* No significant changes.

0.1.31
------

* No significant changes.

0.1.30
------

* ``cpa`` command has been deprecated.

0.1.29
------

* Added ``unicov`` tool.

0.1.28
------

* ``plotcov`` command has been deprecated.

0.1.27
------

* Bumped up the minimally required version of ``vcfgo`` (0.0.3 to 0.0.4).
* Added unit tests.

0.1.26
------

* Integrated Python module ``vcfgo`` to ``peek`` command.
* Integrated Python module ``vcfgo`` to ``viewsnp`` command.

0.1.25
------

* Integrated Python module ``vcfgo`` to ``minivcf`` command.
* Integrated Python module ``vcfgo`` to ``merge`` command.
* Renamed merge.py to mergevcf.py.

0.1.24
------

* Integrated Python module ``vcfgo`` to ``compvcf`` command.

0.1.23
------

* [#4] Fixed bug when using ``bam2gdf`` for the GSTT1 gene with hg38 data.

0.1.22
------

* No significant changes.

0.1.21
------

* No significant changes.

0.1.20
------

* No significant changes.

0.1.19
------

* No significant changes.

0.1.18
------

* No significant changes.

0.1.17
------

* No significant changes.

0.1.16
------

* Updated compvcf.py to also output sample names.

0.1.15
------

* Added ``compvcf`` tool.

0.1.14
------

* No significant changes.

0.1.13
------

* No significant changes.

0.1.12
------

Updated ``compare2`` to also output the target gene column.

0.1.11
------

* Updated ``bam2gt2`` to accept a list of reference samples.

0.1.10
------

* Added ``compare2`` tool.

0.1.9
-----

* Added ``-q 1`` argument to ``bcftools`` caller in bam2vcf.py.

0.1.8
-----

* No significant changes.

0.1.7
-----

* No significant changes.

0.1.6
-----

* Updated sgep.py to support multiple target genes. Consequently, xgep.py was removed.
* Renamed sgep.py to bam2gt2.py.

0.1.5
-----

* Updated sgep.py to support both GATK (bam2vcf2.py) and BCFtools (bam2vcf.py).

0.1.4
-----

* Added bam2vcf2.py (a SGE-based version of bam2vcf.py specifically for ``gatk`` caller).

0.1.3
-----

* Added GATK argument ``java_options`` to bam2vcf.py.

0.1.2
-----

* No significant changes.

0.1.1
-----

* Renamed report.py to gt2html.py.
* Renamed remap.py to bam2bam.py.
* Renamed sges.py to bam2html.py

0.1.0
-----

* Renamed genotype.py to bam2gt.py.
* The new bam2vcf.py can support both GATK and BCFtools for SNP calling.
* The new sgep.py and xgep.py can support both GATK and BCFtools for SNP calling. Consequently, sgea.py and xgea.py were removed.

0.0.80
------

* Added bam2vcf2.py (which will replace bam2vcf.py in near future).

0.0.79
------

* Updated sgep.py, xgep.py and sges.py.

0.0.78
------

* Replaced HaplotypeCaller with BCFtools for sges.py and sgep.py.

0.0.77
------

* Added ``xgea`` tool.

0.0.76
------

* Added ``xgep`` tool.

0.0.75
------

Fixed a bug in the Drugs section for report.py.

0.0.74
------

* Fixed incorrect argument setting for BAM files.

0.0.73
------

* Updated report.py and gt2pt.py.

0.0.72
------

* Updated summary.py and meta.py.

0.0.71
------

* Updated gt2pt.py for CYP2C19 gene.

0.0.70
------

* Added ``gt2pt`` tool (only supports CYP2D6 gene for now).

0.0.69
------

* Updated bam2gdf.py and bam2vcf.py.

0.0.68
------

* Updated sgea.py.

0.0.67
------

* Updated sgep.py.

0.0.66
------

* Updated sges.py and report.py.

0.0.65
------

* Updated genotype.py.

0.0.64
------

* Updated genotype.py.

0.0.63
------

* Added ``genotype`` tool.

0.0.62
------

* Updated bam2vcf.py.

0.0.61
------

* Updated bam2vcf.py.

0.0.60
------

* Added elapsed run time to logging.

0.0.59
------

* Added ``bam2vcf`` tool.

0.0.58
------

* No significant changes.

0.0.57
------

* Updated report.py.

0.0.56
------

* Updated fq2bam.py and remap.py.

0.0.55
------

* No significant changes.

0.0.54
------

* Increased compatibility with Stargazer.

0.0.53
------

* Updated sglib.py.

0.0.52
------

* No significant changes.

0.0.51
------

* No significant changes.

0.0.50
------

* No significant changes.

0.0.49
------

* No significant changes.

0.0.48
------

* Updated ``bam2gdf`` tool to support hg38.

0.0.47
------

* Updated configuration parameters.

0.0.46
------

* Added VCF only mode to ``sges`` tool.

0.0.45
------

* Added VCF only mode to ``sgea`` tool.

0.0.44
------

* Added VCF only mode to ``sgep`` tool.

0.0.43
------

* No significant changes.

0.0.42
------

* Added sglib.py.

0.0.41
------

* No significant changes.

0.0.40
------

* No significant changes.

0.0.39
------

* Added ``snp`` tool.

0.0.38
------

* Added ``peek`` tool.

0.0.37
------

* Added ``liftover`` tool.

0.0.36
------

* Added ``check`` tool.

0.0.35
------

* Added ``plotcov`` tool.

0.0.34
------

* No significant changes.

0.0.33
------

* Added ``cpa`` tool.

0.0.32
------

* Added ``sges`` tool.

0.0.31
------

* Added ``sgep`` tool.

0.0.30
------

* Added ``sgea`` tool.

0.0.29
------

* Added ``fq2bam`` tool.

0.0.28
------

* Added ``remap`` tool.

0.0.27
------

* Added ``compare`` tool.

0.0.26
------

* No significant changes.

0.0.25
------

* Added ``meta`` tool.

0.0.24
------

* Added ``summary`` tool.

0.0.23
------

* No significant changes.

0.0.22
------

* No significant changes.

0.0.21
------

* No significant changes.

0.0.20
------

* Added version.py.

0.0.19
------

* Updated ``VCFFile`` class.

0.0.18
------

* Added ``merge`` tool.

0.0.17
------

* Added ``minivcf`` tool.

0.0.16
------

* No significant changes.

0.0.15
------

* Added Read the Docs.

0.0.14
------

* Added type hints.

0.0.13
------

* Added ``bam2gdf`` tool.

0.0.12
------

* Added ``bam2sdf`` tool.

0.0.11
------

* Added ``sdf2gdf`` tool.

0.0.10
------

* Updated ``pgkb`` tool to be run within Python.

0.0.9
-----

* No significant changes.

0.0.8
-----

* No significant changes.

0.0.7
-----

* Added ``report`` tool.
* Added ``resources`` directory.

0.0.6
-----

* No significant changes.

0.0.5
-----

* No significant changes.

0.0.4
-----

* Added ``pgkb`` tool.

0.0.3
-----

* Added common.py.

0.0.2
-----

* No significant changes.

0.0.1
-----

* Initial release.
