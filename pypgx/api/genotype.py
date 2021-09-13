import warnings

from .. import sdk

import pandas as pd

class GSTT1Genotyper:
    """
    Genotyper for GSTT1.
    """

    def one_row(self, r):
        if r.CNV == 'DeletionHet':
            result = '*A/*0'
        elif r.CNV == 'DeletionHom':
            result = '*0/*0'
        else:
            result = '*A/*A'
        return result

    def genotype(self, df):
        return df.apply(self.one_row, axis=1)

    def __init__(self, df):
        self.results = self.genotype(df)

class UGT2B17Genotyper:
    """
    Genotyper for UGT2B17.
    """

    def one_row(self, r):
        if r.CNV == 'DeletionHet':
            result = '*1/*2'
        elif r.CNV == 'DeletionHom':
            result = '*2/*2'
        else:
            result = '*1/*1'
        return result

    def genotype(self, df):
        return df.apply(self.one_row, axis=1)

    def __init__(self, df):
        self.results = self.genotype(df)

def call_genotypes(alleles=None, cnv_calls=None):
    """
    Call genotypes for target gene.

    Parameters
    ----------
    alleles : pypgx.Archive, optional
        Archive file with the semantic type SampleTable[Alleles].
    cnv_calls : pypgx.Archive, optional
        Archive file with the semantic type SampleTable[CNVCalls].

    Returns
    -------
    pypgx.Archive
        Archive file with the semantic type SampleTable[Genotypes].
    """
    genotypers = {
        'GSTT1': GSTT1Genotyper,
        'UGT2B17': UGT2B17Genotyper,
    }

    # Check the input files.
    if isinstance(alleles, str):
        alleles = sdk.Archive.from_file(alleles)

    if alleles is not None:
        alleles.check('SampleTable[Alleles]')
        def one_row(r):
            r.Haplotype1 = r.Haplotype1.strip(';').split(';')
            r.Haplotype2 = r.Haplotype2.strip(';').split(';')
            return r
        alleles.data = alleles.data.apply(one_row, axis=1)

    if isinstance(cnv_calls, str):
        cnv_calls = sdk.Archive.from_file(cnv_calls)

    if cnv_calls is not None:
        cnv_calls.check('SampleTable[CNVCalls]')

    if alleles is not None and cnv_calls is not None:
        if set(alleles.data.index) != set(cnv_calls.data.index):
            raise ValueError('SampleTable[Alleles] and '
                'SampleTable[CNVCalls] have different samples')
        if alleles.metadata['Gene'] != cnv_calls.metadata['Gene']:
            raise ValueError('Found two different target genes')
        df = pd.concat([alleles.data, cnv_calls.data], axis=1)
        gene = alleles.metadata['Gene']
        assembly = alleles.metadata['Assembly']
    elif alleles is not None and cnv_calls is None:
        df = alleles.data.copy()
        gene = alleles.metadata['Gene']
        assembly = alleles.metadata['Assembly']
    elif alleles is None and cnv_calls is not None:
        df = cnv_calls.data.copy()
        gene = cnv_calls.metadata['Gene']
        assembly = cnv_calls.metadata['Assembly']
    else:
        raise ValueError('Either SampleTable[Alleles] or '
            'SampleTable[CNVCalls] must be provided')

    df = genotypers[gene](df).results.to_frame()
    df.columns = ['Genotype']

    metadata = {}
    metadata['Gene'] = gene
    metadata['Assembly'] = assembly
    metadata['SemanticType'] = 'SampleTable[Genotypes]'

    return sdk.Archive(metadata, df)
