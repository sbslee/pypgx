import os
import io
import zipfile
import tempfile
import copy
import pickle

import pandas as pd
from fuc import pyvcf, pycov, common

class SemanticTypeNotFoundError(Exception):
    """Raise when specified semantic type is not supported."""

class IncorrectSemanticTypeError(Exception):
    """Raise when specified semantic type is incorrect."""

class Archive:
    """
    Class for storing various data.

    Parameters
    ----------
    metadata : dict
        List of metadata lines.
    data : data, results, or model
        Data, results, or model.
    """

    def __init__(self, metadata, data):
        self.metadata = metadata
        self.data = data

    def copy_metadata(self):
        """dict : Copy of the metadata."""
        return copy.deepcopy(self.metadata)

    def to_file(self, fn):
        """
        Create a ZIP file for the Archive.

        Parameters
        ----------
        fn : str
            ZIP file.
        """
        with tempfile.TemporaryDirectory() as t:
            with open(f'{t}/metadata.txt', 'w') as f:
                for k, v in self.metadata.items():
                    if k == 'SemanticType':
                        semantic_type = v
                    f.write(f'{k}={v}\n')
            if 'CovFrame' in self.metadata['SemanticType']:
                self.data.to_file(f'{t}/data.tsv')
            elif 'SampleTable' in self.metadata['SemanticType']:
                self.data.to_csv(f'{t}/data.tsv', sep='\t')
            elif 'VcfFrame' in self.metadata['SemanticType']:
                self.data.to_file(f'{t}/data.vcf')
            elif 'Model' in self.metadata['SemanticType']:
                pickle.dump(self.data, open(f'{t}/data.sav', 'wb'))
            else:
                raise SemanticTypeNotFoundError(self.metadata['SemanticType'])
            zipf = zipfile.ZipFile(fn, 'w', zipfile.ZIP_DEFLATED)
            for root, dirs, files in os.walk(t):
                for file in files:
                    if file == '.DS_Store':
                        continue
                    zipf.write(os.path.join(root, file),
                               os.path.relpath(os.path.join(root, file),
                                               os.path.join(t, '..')))
            zipf.close()

            common.color_print(f'Saved {semantic_type} to: {fn}')

    @classmethod
    def from_file(cls, fn):
        """
        Construct Archive from a ZIP file.

        Parameters
        ----------
        fn : str
            ZIP file.
        """
        metadata = {}
        zf = zipfile.ZipFile(fn)
        parent = zf.filelist[0].filename.split('/')[0]
        with zf.open(f'{parent}/metadata.txt') as f:
            for line in f:
                fields = line.decode('utf-8').strip().split('=')
                metadata[fields[0]] = fields[1]
        if 'CovFrame' in metadata['SemanticType']:
            with zf.open(f'{parent}/data.tsv') as fh:
                data = pycov.CovFrame.from_file(fh)
        elif 'SampleTable' in metadata['SemanticType']:
            with zf.open(f'{parent}/data.tsv') as fh:
                data = pd.read_table(fh, index_col=0)
        elif 'VcfFrame' in metadata['SemanticType']:
            with zf.open(f'{parent}/data.vcf') as fh:
                data = pyvcf.VcfFrame.from_file(fh)
        elif 'Model' in metadata['SemanticType']:
            with zf.open(f'{parent}/data.sav') as fh:
                data = pickle.load(fh)
        else:
            raise SemanticTypeNotFoundError(metadata['SemanticType'])
        return cls(metadata, data)

    def check(self, semantic):
        if self.metadata['SemanticType'] != semantic:
            raise IncorrectSemanticTypeError(f"Expected {semantic}, but found {self.metadata['SemanticType']}")

def zipdir(dir, output):
    """
    Create a ZIP archive of a directory.

    Parameters
    ----------
    dir : str
        Input directory.
    output : str
        Output file.
    """
    zipf = zipfile.ZipFile(output, 'w', zipfile.ZIP_DEFLATED)
    for root, dirs, files in os.walk(dir):
        for file in files:
            if file == '.DS_Store':
                continue
            zipf.write(os.path.join(root, file),
                       os.path.relpath(os.path.join(root, file),
                                       os.path.join(dir, '..')))
    zipf.close()

def parse_pharmvar(fn):
    """
    Parse PharmVar gene data.

    Parameters
    ----------
    fn : str
        Gene data directory.
    """
    gene = os.path.basename(fn).split('-')[0]

    rs_dict = {}

    vfs = {'GRCh37': [], 'GRCh38': []}
    alleles = {}
    for i, assembly in enumerate(['GRCh37', 'GRCh38']):
        for r, d, f in os.walk(f'{fn}/{assembly}'):
            for file in f:
                if file.endswith('.tsv'):
                    df = pd.read_table(f'{r}/{file}', comment='#')
                if file.endswith('.vcf'):
                    vf = pyvcf.VcfFrame.from_file(f'{r}/{file}')
                    vfs[assembly].append(vf)
        chrom = vfs['GRCh37'][0].contigs[0]
        for j, r in df.iterrows():
            name = r['Haplotype Name'].replace(gene, '')
            if name not in alleles:
                alleles[name] = [[], []]
            if pd.isna(r['Variant Allele']):
                continue
            variant = f"{chrom}-{r['Variant Start']}-{r['Reference Allele']}-{r['Variant Allele']}"
            rs_dict[variant] = r['rsID']
            alleles[name][i].append(variant)

    variants = {'GRCh37': {}, 'GRCh38': {}}

    for name in alleles:
        for i, assembly in enumerate(['GRCh37', 'GRCh38']):
            for variant in alleles[name][i]:
                if variant not in variants[assembly]:
                    variants[assembly][variant] = []
                if name not in variants[assembly][variant]:
                    variants[assembly][variant].append(name)

    for name in alleles:
        alleles[name] = [','.join(alleles[name][0]), ','.join(alleles[name][1])]

    df1 = pd.DataFrame(alleles).T
    df1.columns = ['GRCh37', 'GRCh38']
    df1 = df1.replace('', 'N/A')
    df1.to_csv(f'{gene}-allele-table.csv')

    def func(r):
        if len(r.REF) == len(r.ALT) == 1:
            return f'{r.CHROM}-{r.POS}-{r.REF}-{r.ALT}'
        elif len(r.REF) == 1 and len(r.ALT) > 1:
            return f'{r.CHROM}-{r.POS}---{r.ALT[1:]}'
        elif len(r.REF) > 1 and len(r.ALT) == 1:
            return f'{r.CHROM}-{r.POS+1}-{r.REF[1:]}--'
        else:
            raise ValueError('Something went wrong')

    for assembly in ['GRCh37', 'GRCh38']:
        df2 = pyvcf.merge(vfs[assembly]).chr_prefix().df
        df2['Name'] = df2.apply(func, axis=1)
        df2['Alleles'] = df2.apply(lambda r: ','.join(variants[assembly][r.Name]), axis=1)
        df2['rsID'] = df2.apply(lambda r: rs_dict[r.Name], axis=1)
        df2.to_csv(f'{gene}-{assembly}.csv')
