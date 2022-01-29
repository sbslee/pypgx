from setuptools import setup, find_packages

exec(open('pypgx/version.py').read())

requirements = ['fuc', 'scikit-learn']

setup(
    name='pypgx',
    version=__version__,
    author='Seung-been "Steven" Lee',
    author_email='sbstevenlee@gmail.com',
    description='A Python package for pharmacogenomics (PGx) research',
    long_description=open('README.rst').read(),
    long_description_content_type='text/x-rst',
    url='https://github.com/sbslee/pypgx',
    packages=find_packages(),
    license='MIT',
    install_requires=requirements,
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License'
    ],
    package_data={
        'pypgx.api': ['beagle.28Jun21.220.jar', 'data/*']
    },
    entry_points={'console_scripts': ['pypgx=pypgx.__main__:main']}
)
