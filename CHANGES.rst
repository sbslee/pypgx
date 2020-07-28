Changelog
*********

v0.1.8
------

* No significant changes.

v0.1.7
------

* No significant changes.

v0.1.6
------

* Updated sgep.py to support multiple target genes. Consequently, xgep.py 
  was removed.
* Renamed sgep.py to bam2gt2.py.

v0.1.5
------

* Updated sgep.py to support both GATK (bam2vcf2.py) and BCFtools 
  (bam2vcf.py).

v0.1.4
------

* Added bam2vcf2.py (a SGE-based version of bam2vcf.py specifically for 
  ``gatk`` caller).

v0.1.3
------

* Added GATK argument ``java_options`` to bam2vcf.py.

v0.1.2
------

* No significant changes.

v0.1.1
------

* Renamed report.py to gt2html.py.
* Renamed remap.py to bam2bam.py.
* Renamed sges.py to bam2html.py

v0.1.0
------

* Renamed genotype.py to bam2gt.py.
* The new bam2vcf.py can support both GATK and BCFtools for SNP calling.
* The new sgep.py and xgep.py can support both GATK and BCFtools for SNP 
  calling. Consequently, sgea.py and xgea.py were removed.

v0.0.80
-------

* Added bam2vcf2.py (which will replace bam2vcf.py in near future).

v0.0.79
-------

* Updated sgep.py, xgep.py and sges.py.

v0.0.78
-------

* Replaced HaplotypeCaller with BCFtools for sges.py and sgep.py.

v0.0.77
-------

* Added ``xgea`` tool.

v0.0.76
-------

* Added ``xgep`` tool.

v0.0.75
-------

Fixed a bug in the Drugs section for report.py.

v0.0.74
-------

* Fixed incorrect argument setting for BAM files.

v0.0.73
-------

* Updated report.py and gt2pt.py.

v0.0.72
-------

* Updated summary.py and meta.py.

v0.0.71
-------

* Updated gt2pt.py for CYP2C19 gene.

v0.0.70
-------

* Added ``gt2pt`` tool (only supports CYP2D6 gene for now).

v0.0.69
-------

* Updated bam2gdf.py and bam2vcf.py.

v0.0.68
-------

* Updated sgea.py.

v0.0.67
-------

* Updated sgep.py.

v0.0.66
-------

* Updated sges.py and report.py.

v0.0.65
-------

* Updated genotype.py.

v0.0.64
-------

* Updated genotype.py.

v0.0.63
-------

* Added ``genotype`` tool.

v0.0.62
-------

* Updated bam2vcf.py.

v0.0.61
-------

* Updated bam2vcf.py.

v0.0.60
-------

* Added elapsed run time to logging.

v0.0.59
-------

* Added ``bam2vcf`` tool.

v0.0.58
-------

* No significant changes.

v0.0.57
-------

* Updated report.py.

v0.0.56
-------

* Updated fq2bam.py and remap.py.

v0.0.55
-------

* No significant changes.

v0.0.54
-------

* Increased compatibility with Stargazer.

v0.0.53
-------

* Updated sglib.py.

v0.0.52
-------

* No significant changes.

v0.0.51
-------

* No significant changes.

v0.0.50
-------

* No significant changes.

v0.0.49
-------

* No significant changes.

v0.0.48
-------

* Updated ``bam2gdf`` tool to support hg38.

v0.0.47
-------

* Updated configuration parameters.

v0.0.46
-------

* Added VCF only mode to ``sges`` tool.

v0.0.45
-------

* Added VCF only mode to ``sgea`` tool.

v0.0.44
-------

* Added VCF only mode to ``sgep`` tool.

v0.0.43
-------

* No significant changes.

v0.0.42
-------

* Added sglib.py.

v0.0.41
-------

* No significant changes.

v0.0.40
-------

* No significant changes.

v0.0.39
-------

* Added ``snp`` tool.

v0.0.38
-------

* Added ``peek`` tool.

v0.0.37
-------

* Added ``liftover`` tool.

v0.0.36
-------

* Added ``check`` tool.

v0.0.35
-------

* Added ``plotcov`` tool.

v0.0.34
-------

* No significant changes.

v0.0.33
-------

* Added ``cpa`` tool.

v0.0.32
-------

* Added ``sges`` tool.


v0.0.31
-------

* Added ``sgep`` tool.

v0.0.30
-------

* Added ``sgea`` tool.

v0.0.29
-------

* Added ``fq2bam`` tool.

v0.0.28
-------

* Added ``remap`` tool.

v0.0.27
-------

* Added ``compare`` tool.

v0.0.26
-------

* No significant changes.

v0.0.25
-------

* Added ``meta`` tool.

v0.0.24
-------

* Added ``summary`` tool.

v0.0.23
-------

* No significant changes.

v0.0.22
-------

* No significant changes.

v0.0.21
-------

* No significant changes.

v0.0.20
-------

* Added version.py.

v0.0.19
-------

* Updated ``VCFFile`` class.

v0.0.18
-------

* Added ``merge`` tool.

v0.0.17
-------

* Added ``minivcf`` tool.

v0.0.16
-------

* No significant changes.

v0.0.15
-------

* Added Read the Docs.

v0.0.14
-------

* Added type hints.

v0.0.13
-------

* Added ``bam2gdf`` tool.

v0.0.12
-------

* Added ``bam2sdf`` tool.

v0.0.11
-------

* Added ``sdf2gdf`` tool.

v0.0.10
-------

* Updated ``pgkb`` tool to be run within Python.

v0.0.9
------

* No significant changes.

v0.0.8
------

* No significant changes.

v0.0.7
------

* Added ``report`` tool.
* Added ``resources`` directory.

v0.0.6
------

* No significant changes.

v0.0.5
------

* No significant changes.

v0.0.4
------

* Added ``pgkb`` tool.

v0.0.3
------

* Added common.py.

v0.0.2
------

* No significant changes.

v0.0.1
------

* Initial release.
