import warnings

from .. import sdk

import pandas as pd

class UGT2B17Genotyper:
    """
    Genotyper for UGT2B17.
    """

    def one_row(self, r):
        if r.CNV == 'DeletionHet':
            r['Genotype'] = '*1/*2'
        elif r.CNV == 'DeletionHom':
            r['Genotype'] = '*2/*2'
        else:
            r['Genotype'] = '*1/*1'
        return r

    def genotype(self, df):
        return df.apply(self.one_row, axis=1)

    def __init__(self, df):
        self.results = self.genotype(df)

def call_genotypes(alleles=None, cnv=None):
    """
    Call genotypes for specified gene.

    Parameters
    ----------
    alleles : pypgx.Archive, optional
        Archive file with the semantic type SampleTable[Alleles].
    cnv : pypgx.Archive, optional
        Archive file with the semantic type SampleTable[CNVCalls].

    Returns
    -------
    pypgx.Archive
        Archive file with the semantic type SampleTable[Genotypes].
    """
    genotypers = {
        'UGT2B17': UGT2B17Genotyper,
    }

    # Check the inputs.
    if alleles is not None and isinstance(alleles, str):
        alleles = sdk.Archive.from_file(alleles)

    if alleles is not None and isinstance(cnv, str):
        cnv = sdk.Archive.from_file(cnv)

    alleles.check('SampleTable[Alleles]')

    cnv.check('SampleTable[CNVCalls]')

    if set(alleles.data.index) != set(cnv.data.index):
        raise ValueError('SampleTable[Alleles] and SampleTable[CNVCalls] have different samples')

    def one_row(r):
        r.Haplotype1 = r.Haplotype1.strip(';').split(';')
        r.Haplotype2 = r.Haplotype2.strip(';').split(';')
        return r
    alleles.data = alleles.data.apply(one_row, axis=1)

    if alleles.metadata['Gene'] != cnv.metadata['Gene']:
        raise ValueError('Found two different target genes')

    df = pd.concat([alleles.data, cnv.data], axis=1)

    if df.isna().any().any():
        m = (f'Will drop {df.isna().any(axis=1).sum()} samples because of '
            'missing data')
        df = df.dropna()
        warnings.warn(m)

    df = genotypers[alleles.metadata['Gene']](df).results
    def one_row(r):
        r.Haplotype1 = ';'.join(r.Haplotype1) + ';'
        r.Haplotype2 = ';'.join(r.Haplotype2) + ';'
        return r
    df = df.apply(one_row, axis=1)

    metadata = {}
    metadata['Gene'] = alleles.metadata['Gene']
    metadata['Assembly'] = alleles.metadata['Assembly']
    metadata['SemanticType'] = 'SampleTable[Genotypes]'

    return sdk.Archive(metadata, df)
