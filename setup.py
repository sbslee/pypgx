from setuptools import setup, find_packages

with open("README.md", "r") as readme_file:
    readme = readme_file.read()

requirements = ["requests>=2"]

setup(
    name="pypgx",
    version="0.0.2",
    author='Seung-been "Steven" Lee',
    author_email="sbstevenlee@gmail.com",
    description="A Python package for PGx research",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/pypgx/pypgx/",
    packages=find_packages(),
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ],
    entry_points={"console_scripts": ["pypgx=pypgx.__main__:main"]}
)
