import pydoc
import subprocess

import pypgx
import pypgx.api
from pypgx.api import core
from pypgx.cli import commands

submodules = ['core', 'genotype', 'pipeline', 'plot', 'utils']

credit = """
..
   This file was automatically generated by docs/create.py.
"""

pypgx_help = subprocess.run(['pypgx', '-h'], capture_output=True, text=True, check=True).stdout
pypgx_help = '\n'.join(['   ' + x for x in pypgx_help.splitlines()])

submodule_help = ''

for submodule in submodules:
    description = pydoc.getdoc(getattr(pypgx.api, submodule)).split('\n\n')[0].replace('\n', ' ')
    submodule_help += f'- **{submodule}** : {description}\n'

d = dict(credit=credit, pypgx_help=pypgx_help, submodule_help=submodule_help)

# -- README.rst ---------------------------------------------------------------

readme = """
{credit}
README
******

.. image:: https://badge.fury.io/py/pypgx.svg
    :target: https://badge.fury.io/py/pypgx

.. image:: https://readthedocs.org/projects/pypgx/badge/?version=latest
    :target: https://pypgx.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://anaconda.org/bioconda/pypgx/badges/version.svg
   :target: https://anaconda.org/bioconda/pypgx

.. image:: https://anaconda.org/bioconda/pypgx/badges/license.svg
   :target: https://github.com/sbslee/pypgx/blob/master/LICENSE

.. image:: https://anaconda.org/bioconda/pypgx/badges/downloads.svg
   :target: https://anaconda.org/bioconda/pypgx/files

.. image:: https://anaconda.org/bioconda/pypgx/badges/installer/conda.svg
   :target: https://conda.anaconda.org/bioconda

Introduction
============

The main purpose of the PyPGx package is to provide a unified platform for pharmacogenomics (PGx) research.

The package is written in Python, and supports both command line interface (CLI) and application programming interface (API) whose documentations are available at the `Read the Docs <https://pypgx.readthedocs.io/en/latest/>`_.

Your contributions (e.g. feature ideas, pull requests) are most welcome.

| Author: Seung-been "Steven" Lee
| Email: sbstevenlee@gmail.com
| License: MIT License

Installation
============

Following packages are required to run PyPGx:

.. list-table::
   :header-rows: 1

   * - Package
     - Anaconda
     - PyPI
   * - ``fuc``
     - ✅
     - ✅
   * - ``scikit-learn``
     - ✅
     - ✅
   * - ``openjdk``
     - ✅
     - ❌

There are various ways you can install PyPGx. The recommended way is via conda (`Anaconda <https://www.anaconda.com/>`__):

.. code-block:: text

   $ conda install -c bioconda pypgx

Above will automatically download and install all the dependencies as well. Alternatively, you can use pip (`PyPI <https://pypi.org/>`__) to install PyPGx and all of its dependencies except ``openjdk`` (i.e. Java JDK must be installed separately):

.. code-block:: text

   $ pip install pypgx

Finally, you can clone the GitHub repository and then install PyPGx locally:

.. code-block:: text

   $ git clone https://github.com/sbslee/pypgx
   $ cd pypgx
   $ pip install .

The nice thing about this approach is that you will have access to development versions that are not available in Anaconda or PyPI. For example, you can access a development branch with the ``git checkout`` command. When you do this, please make sure your environment already has all the dependencies installed.

Archive file, semantic type, and metadata
=========================================

In order to efficiently store and transfer data, PyPGx uses the ZIP archive file format (``.zip``) which supports lossless data compression. Each archive file created by PyPGx has a metadata file (``metadata.txt``) and a data file (e.g. ``data.tsv``, ``data.vcf``). A metadata file contains important information about the data file within the same archive, which is expressed as pairs of ``=``-separated keys and values (e.g. ``Assembly=GRCh37``):

.. list-table::
    :widths: 20 40 40
    :header-rows: 1

    * - Metadata
      - Description
      - Examples
    * - ``Assembly``
      - Reference genome assembly.
      - ``GRCh37``, ``GRCh38``
    * - ``Control``
      - Control gene.
      - ``VDR``, ``chr1:10000-20000``
    * - ``Gene``
      - Target gene.
      - ``CYP2D6``, ``GSTT1``
    * - ``Platform``
      - Genotyping platform.
      - ``WGS``, ``Targeted``, ``Chip``
    * - ``Program``
      - Name of the phasing program.
      - ``Beagle``
    * - ``Samples``
      - Samples used for inter-sample normalization.
      - ``NA07000,NA10854,NA11993``
    * - ``SemanticType``
      - Semantic type of the archive.
      - ``CovFrame[CopyNumber]``, ``Model[CNV]``

Notably, all archive files have defined semantic types, which allows us to ensure that the data that is passed to a PyPGx command (CLI) or method (API) is meaningful for the operation that will be performed. Below is a list of currently defined semantic types:

- ``CovFrame[CopyNumber]``
    * CovFrame for storing target gene's per-base copy number which is computed from read depth with control statistics.
    * Requires following metadata: ``Gene``, ``Assembly``, ``SemanticType``, ``Platform``, ``Control``, ``Samples``.
- ``CovFrame[DepthOfCoverage]``
    * CovFrame for storing read depth for all target genes with SV.
    * Requires following metadata: ``Assembly``, ``SemanticType``, ``Platform``.
- ``CovFrame[ReadDepth]``
    * CovFrame for storing read depth for single target gene.
    * Requires following metadata: ``Gene``, ``Assembly``, ``SemanticType``, ``Platform``.
- ``Model[CNV]``
    * Model for calling CNV in target gene.
    * Requires following metadata: ``Gene``, ``Assembly``, ``SemanticType``, ``Control``.
- ``SampleTable[Alleles]``
    * TSV file for storing target gene's candidate star alleles for each sample.
    * Requires following metadata: ``Platform``, ``Gene``, ``Assembly``, ``SemanticType``, ``Program``.
- ``SampleTable[CNVCalls]``
    * TSV file for storing target gene's CNV call for each sample.
    * Requires following metadata: ``Gene``, ``Assembly``, ``SemanticType``, ``Control``.
- ``SampleTable[Genotypes]``
    * TSV file for storing target gene's genotype call for each sample.
    * Requires following metadata: ``Gene``, ``Assembly``, ``SemanticType``.
- ``SampleTable[Phenotypes]``
    * TSV file for storing target gene's phenotype call for each sample.
    * Requires following metadata: ``Gene``, ``SemanticType``.
- ``SampleTable[Results]``
    * TSV file for storing various results for each sample.
    * Requires following metadata: ``Gene``, ``Assembly``, ``SemanticType``.
- ``SampleTable[Statistcs]``
    * TSV file for storing control gene's various statistics on read depth for each sample. Used for converting target gene's read depth to copy number.
    * Requires following metadata: ``Control``, ``Assembly``, ``SemanticType``, ``Platform``.
- ``VcfFrame[Consolidated]``
    * VcfFrame for storing target gene's consolidated variant data.
    * Requires following metadata: ``Platform``, ``Gene``, ``Assembly``, ``SemanticType``, ``Program``.
- ``VcfFrame[Imported]``
    * VcfFrame for storing target gene's raw variant data.
    * Requires following metadata: ``Platform``, ``Gene``, ``Assembly``, ``SemanticType``.
- ``VcfFrame[Phased]``
    * VcfFrame for storing target gene's phased variant data.
    * Requires following metadata: ``Platform``, ``Gene``, ``Assembly``, ``SemanticType``, ``Program``.

Getting help
============
For detailed documentations on the CLI and API, please refer to the `Read the Docs <https://pypgx.readthedocs.io/en/latest/>`_.

For getting help on the CLI:

.. code-block:: text

   $ pypgx -h

{pypgx_help}

For getting help on a specific command (e.g. call-genotypes):

.. code-block:: text

   $ pypgx call-genotypes -h

Below is the list of submodules available in the API:

{submodule_help}

For getting help on a specific submodule (e.g. utils):

.. code:: python3

   >>> from pypgx.api import utils
   >>> help(utils)

CLI examples
============

We can print the metadata of an archive file:

.. code-block:: text

    $ pypgx print-metadata CYP2D6-copy-number.zip

Above will print:

.. code-block:: text

    Gene=CYP2D6
    Assembly=GRCh37
    SemanticType=CovFrame[CopyNumber]
    Platform=WGS
    Control=VDR
    Samples=None

We can run the NGS pipeline for the *CYP2D6* gene:

.. code-block:: text

    $ pypgx run-ngs-pipeline \\
    CYP2D6 \\
    CYP2D6-pipeline \\
    --variants variants.vcf \\
    --depth-of-coverage depth-of-coverage.zip \\
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

API examples
============

We can obtain allele function for the *CYP2D6* gene:

.. code:: python3

    >>> import pypgx
    >>> pypgx.get_function('CYP2D6', '*1')
    'Normal Function'
    >>> pypgx.get_function('CYP2D6', '*4')
    'No Function'
    >>> pypgx.get_function('CYP2D6', '*22')
    'Uncertain Function'
    >>> pypgx.get_function('CYP2D6', '*140')
    'Unknown Function'

We can predict phenotype for the *CYP2D6* gene based on two haplotype calls:

.. code:: python3

    >>> import pypgx
    >>> pypgx.predict_phenotype('CYP2D6', '*4', '*5')   # Both alleles have no function
    'Poor Metabolizer'
    >>> pypgx.predict_phenotype('CYP2D6', '*5', '*4')   # The order of alleles does not matter
    'Poor Metabolizer'
    >>> pypgx.predict_phenotype('CYP2D6', '*1', '*22')  # *22 has uncertain function
    'Indeterminate'
    >>> pypgx.predict_phenotype('CYP2D6', '*1', '*1x2') # Gene duplication
    'Ultrarapid Metabolizer'
""".format(**d)

readme_file = f'{core.PROGRAM_PATH}/README.rst'

with open(readme_file, 'w') as f:
    f.write(readme.lstrip())

# -- cli.rst -----------------------------------------------------------------

cli = """
{credit}

CLI
***

Introduction
============

This section describes the command line interface (CLI) for PyPGx.

For getting help on the CLI:

.. code-block:: text

   $ pypgx -h

{pypgx_help}

For getting help on a specific command (e.g. call-genotypes):

.. code-block:: text

   $ pypgx call-genotypes -h

""".format(**d)

for command in commands:
    s = f'{command}\n'
    s += '=' * (len(s)-1) + '\n'
    s += '\n'
    s += '.. code-block:: text\n'
    s += '\n'
    s += f'   $ pypgx {command} -h\n'
    command_help = subprocess.run(['pypgx', command, '-h'], capture_output=True, text=True, check=True).stdout
    command_help = '\n'.join(['   ' + x for x in command_help.splitlines()])
    s += command_help + '\n'
    s += '\n'
    cli += s

cli_file = f'{core.PROGRAM_PATH}/docs/cli.rst'

with open(cli_file, 'w') as f:
    f.write(cli.lstrip())

# -- api.rst -----------------------------------------------------------------

api = """
{credit}
API
***

Introduction
============

This section describes the application programming interface (API) for PyPGx.

Below is the list of submodules available in the API:

{submodule_help}

For getting help on a specific submodule (e.g. utils):

.. code:: python3

   from pypgx.api import utils
   help(utils)

""".format(**d)

for submodule in submodules:
    s = f'{submodule}\n'
    s += '=' * (len(s)-1) + '\n'
    s += '\n'
    s += f'.. automodule:: pypgx.api.{submodule}\n'
    s += '   :members:\n'
    s += '\n'
    api += s

with open(f'{core.PROGRAM_PATH}/docs/api.rst', 'w') as f:
    f.write(api.lstrip())

# -- sdk.rst -----------------------------------------------------------------

sdk = """
{credit}
SDK
***

Introduction
============

This section describes the software development kit (SDK) for PyPGx.

utils
=====

.. automodule:: pypgx.sdk.utils
   :members:

""".format(**d)

with open(f'{core.PROGRAM_PATH}/docs/sdk.rst', 'w') as f:
    f.write(sdk.lstrip())
