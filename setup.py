from setuptools import setup, find_packages

with open("README.rst", "r") as readme_file:
    readme = readme_file.read()

exec(open("pypgx/version.py").read())

requirements = ["requests>=2", "pandas>=1.0.0", "bs4>=0.0.1",
    "lxml>=4.5.0", "pysam>=0.16.0",]

setup(
    name="pypgx",
    version=__version__,
    author='Seung-been "Steven" Lee',
    author_email="sbstevenlee@gmail.com",
    description=(
        "PyPGx is a Python package for pharmacogenomics (PGx) research"),
    long_description=readme,
    long_description_content_type="text/x-rst",
    url="https://github.com/pypgx/pypgx/",
    packages=find_packages(),
    package_data={
        "pypgx.resources.cpic": ["cpicPairs.csv"],
        "pypgx.resources.sg": ["gene_table.txt"],
        "pypgx.resources.pgkb": ["actions.txt"],
    },
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ],
    entry_points={"console_scripts": ["pypgx=pypgx.__main__:main"]}
)
