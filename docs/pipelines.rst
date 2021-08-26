Pipelines
*********

This section describes various pipelines for PyPGx.

Statistical phasing
===================

Below is a non-exhaustive list of popular and publicly available tools used for statistical phasing:

- `Beagle <https://faculty.washington.edu/browning/beagle/beagle.html>`__
- `SHAPEIT <https://odelaneau.github.io/shapeit4/>`__
- `Minimac <https://github.com/statgen/Minimac4>`__
- `Eagle <https://alkesgroup.broadinstitute.org/Eagle/>`__

For Beagle, you can download a specific version from the website (e.g. 28Jun21.220):

.. code-block:: text

    $ wget https://faculty.washington.edu/browning/beagle/beagle.28Jun21.220.jar

Each tool estimates haplotype phase either within a genotyped cohort or using a phased reference panel.

If your input data is GRCh37, I recommend using the 1000 Genomes Project phase 3 reference panel. You can easily download it from the link provided by the authors of Beagle:

.. code-block:: text

    $ wget -r --no-parent http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/

For running Beagle:

.. code-block:: text

    $ java -Xmx2g -jar beagle.28Jun21.220.jar gt=in.vcf chrom=chr22:42512500-42551883 ref=chr22.1kg.phase3.v5a.vcf.gz out=out.vcf impute=false
