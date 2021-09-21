"""
The genotype submodule is a suite of tools for accurately predicting genotype
calls.
"""

import warnings

from .. import sdk
from . import utils

import pandas as pd

###################
# Private methods #
###################

def _call_duplication(r):
    """
    Call whole gene duplication.
    """
    a1 = r.Haplotype1[0]
    a2 = r.Haplotype2[0]

    if r.VariantData[a1]:
        h1 = all([x > 0.5 for x in r.VariantData[a1][1]])
    else:
        h1 = True

    if r.VariantData[a2]:
        h2 = all([x > 0.5 for x in r.VariantData[a2][1]])
    else:
        h2 = True

    if h1 and h2:
        if a1 == a2:
            result = '/'.join(sorted([a2, a1+'x2']))
        elif r.VariantData[a1] and not r.VariantData[a2]:
            result = '/'.join(sorted([a2, a1+'x2']))
        elif not r.VariantData[a1] and r.VariantData[a2]:
            result = '/'.join(sorted([a1, a2+'x2']))
        else:
            result = 'Unassigned'
    elif h1 and not h2:
        result = '/'.join(sorted([a2, a1+'x2']))
    elif not h1 and h2:
        result = '/'.join(sorted([a1, a2+'x2']))
    else:
        result = 'Unassigned'

    return result

###############################
# Public classes and methods  #
###############################

class SimpleGenotyper:
    """
    Genotyper for genes without SV.
    """

    def one_row(self, r):
        alleles = [r.Haplotype1[0], r.Haplotype2[0]]
        result = '/'.join(sorted(alleles))
        return result

    def __init__(self, df):
        self.results = df.apply(self.one_row, axis=1)

class CYP2A6Genotyper:
    """
    Genotyper for CYP2A6.
    """

    def one_row(self, r):
        alleles = [r.Haplotype1[0], r.Haplotype2[0]]
        if r.CNV == 'Normal':
            result = '/'.join(sorted(alleles))
        elif r.CNV == 'DeletionHom':
            result = '*4/*4'
        elif r.CNV == 'DeletionHet':
            if alleles[0] == alleles[1]:
                result = '/'.join(sorted([alleles[0], '*4']))
            elif alleles[0] == '*1':
                result = '/'.join(sorted([alleles[1], '*4']))
            elif alleles[1] == '*1':
                result = '/'.join(sorted([alleles[0], '*4']))
            else:
                result = 'Unassigned'
        else:
            result = 'Unassigned'
        return result

    def __init__(self, df):
        self.results = df.apply(self.one_row, axis=1)

class CYP2B6Genotyper:
    """
    Genotyper for CYP2B6.
    """

    def one_row(self, r):
        alleles = [r.Haplotype1[0], r.Haplotype2[0]]
        if r.CNV == 'Normal':
            result = '/'.join(sorted(alleles))
        elif r.CNV == 'Hybrid':
            h1 = len(r.Haplotype1)
            h2 = len(r.Haplotype2)
            if h1 == h2:
                result = '/'.join(sorted([alleles[0], '*29']))
            elif h1 < h2:
                result = '/'.join(sorted([alleles[1], '*29']))
            elif h1 > h2:
                result = '/'.join(sorted([alleles[0], '*29']))
            else:
                result = 'Unassigned'
        else:
            result = 'Unassigned'
        return result

    def __init__(self, df):
        self.results = df.apply(self.one_row, axis=1)

class CYP2D6Genotyper:
    """
    Genotyper for CYP2D6.
    """

    def one_row(self, r):
        alleles = [r.Haplotype1[0], r.Haplotype2[0]]
        if r.CNV == 'Normal':
            result = '/'.join(sorted(alleles))
        elif r.CNV == 'DeletionHom':
            result = '*5/*5'
        elif r.CNV == 'DeletionHet':
            if alleles[0] == alleles[1]:
                result = '/'.join(sorted([alleles[0], '*5']))
            elif alleles[0] == '*1':
                result = '/'.join(sorted([alleles[1], '*5']))
            elif alleles[1] == '*1':
                result = '/'.join(sorted([alleles[0], '*5']))
            else:
                result = 'Unassigned'
        elif r.CNV == 'Duplication':
            result = _call_duplication(r)
        elif r.CNV == 'Tandem1':
            h1 = '*4' in r.Haplotype1
            h2 = '*4' in r.Haplotype2
            if h1 and h2:
                result = '/'.join(sorted([alleles[0], '*68+*4']))
            elif h1 and not h2:
                result = '/'.join(sorted([alleles[1], '*68+*4']))
            elif not h1 and h2:
                result = '/'.join(sorted([alleles[0], '*68+*4']))
            else:
                result = 'Unassigned'
        elif r.CNV == 'Tandem2A':
            h1 = '*10' in r.Haplotype1
            h2 = '*10' in r.Haplotype2
            if h1 and h2:
                result = '/'.join(sorted([alleles[0], '*36+*10']))
            elif h1 and not h2:
                result = '/'.join(sorted([alleles[1], '*36+*10']))
            elif not h1 and h2:
                result = '/'.join(sorted([alleles[0], '*36+*10']))
            else:
                result = 'Unassigned'
        elif r.CNV == 'Tandem2B':
            h1 = '*10' in r.Haplotype1
            h2 = '*10' in r.Haplotype2
            if h1 and h2:
                result = '/'.join(sorted(['*36+*10', '*36+*10']))
            elif h1 and not h2:
                result = '/'.join(sorted([alleles[1], '*36x2+*10']))
            elif not h1 and h2:
                result = '/'.join(sorted([alleles[0], '*36x2+*10']))
            else:
                result = 'Unassigned'
        elif 'DeletionHet' in r.CNV and 'Tandem1' in r.CNV:
            if '*4' in r.Haplotype1 or '*4' in r.Haplotype2:
                result = '/'.join(sorted(['*5', '*68+*4']))
            else:
                result = 'Unassigned'
        elif 'Duplication' in r.CNV and 'Tandem1' in r.CNV:
            h1 = '*4' in r.Haplotype1
            h2 = '*4' in r.Haplotype2
            if h1 and h2:
                result = '/'.join(sorted(['*4x2', '*68+*4']))
            elif h1 and not h2:
                result = '/'.join(sorted([alleles[1]+'x2', '*68+*4']))
            elif not h1 and h2:
                result = '/'.join(sorted([alleles[0]+'x2', '*68+*4']))
            else:
                result = 'Unassigned'
        else:
            result = 'Unassigned'

        return result

    def __init__(self, df):
        self.results = df.apply(self.one_row, axis=1)

class CYP2E1Genotyper:
    """
    Genotyper for CYP2E1.
    """

    def one_row(self, r):
        alleles = [r.Haplotype1[0], r.Haplotype2[0]]
        if r.CNV == 'Normal':
            result = '/'.join(sorted(alleles))
        elif r.CNV == 'PartialDuplication':
            h1 = '*4'in r.Haplotype1 and '*7' in r.Haplotype1
            h2 = '*4'in r.Haplotype2 and '*7' in r.Haplotype2
            if h1 and h2:
                result = '/'.join(sorted([alleles[0], '*S1']))
            elif h1 and not h2:
                result = '/'.join(sorted([alleles[1], '*S1']))
            elif not h1 and h2:
                result = '/'.join(sorted([alleles[0], '*S1']))
            else:
                result = 'Unassigned'
        elif r.CNV == 'Duplication':
            h1 = '*7'in r.Haplotype1
            h2 = '*7'in r.Haplotype2
            if h1 and h2:
                result = '/'.join(sorted([alleles[0], '*7x2']))
            elif h1 and not h2:
                result = '/'.join(sorted([alleles[1], '*7x2']))
            elif not h1 and h2:
                result = '/'.join(sorted([alleles[0], '*7x2']))
            else:
                result = 'Unassigned'
        else:
            result = 'Unassigned'
        return result

    def __init__(self, df):
        self.results = df.apply(self.one_row, axis=1)

class GSTM1Genotyper:
    """
    Genotyper for GSTM1.
    """

    def one_row(self, r):
        alleles = [r.Haplotype1[0], r.Haplotype2[0]]
        if r.CNV == 'DeletionHet':
            result = '/'.join(sorted([alleles[0], '*0']))
        elif r.CNV == 'DeletionHom':
            result = '*0/*0'
        elif r.CNV == 'Normal':
            result = '/'.join(sorted(alleles))
        else:
            result = 'Unassigned'
        return result

    def __init__(self, df):
        self.results = df.apply(self.one_row, axis=1)

class GSTT1Genotyper:
    """
    Genotyper for GSTT1.
    """

    def one_row(self, r):
        if r.CNV == 'DeletionHet':
            result = '*A/*0'
        elif r.CNV == 'DeletionHom':
            result = '*0/*0'
        elif r.CNV == 'Normal':
            result = '*A/*A'
        else:
            result = 'Unassigned'
        return result

    def __init__(self, df):
        self.results = df.apply(self.one_row, axis=1)

class SLC22A2Genotyper:
    """
    Genotyper for SLC22A2.
    """

    def one_row(self, r):
        alleles = [r.Haplotype1[0], r.Haplotype2[0]]
        if r.CNV == 'Normal':
            result = '/'.join(sorted(alleles))
        elif r.CNV == 'Intron9Deletion':
            h1 = '*K432Q'in r.Haplotype1
            h2 = '*K432Q'in r.Haplotype2
            if h1 and h2:
                result = '/'.join(sorted([alleles[0], '*S1']))
            elif h1 and not h2:
                result = '/'.join(sorted([alleles[1], '*S1']))
            elif not h1 and h2:
                result = '/'.join(sorted([alleles[0], '*S1']))
            else:
                result = 'Unassigned'
        elif r.CNV == 'Exon11Deletion':
            h1 = r.Haplotype1 == ['*3']
            h2 = r.Haplotype2 == ['*3']
            if h1 and h2:
                result = '/'.join(sorted([alleles[0], '*S2']))
            elif h1 and not h2:
                result = '/'.join(sorted([alleles[1], '*S2']))
            elif not h1 and h2:
                result = '/'.join(sorted([alleles[0], '*S2']))
            else:
                result = 'Unassigned'
        else:
            result = 'Unassigned'
        return result

    def __init__(self, df):
        self.results = df.apply(self.one_row, axis=1)

class UGT2B15Genotyper:
    """
    Genotyper for UGT2B15.
    """

    def one_row(self, r):
        alleles = [r.Haplotype1[0], r.Haplotype2[0]]
        if r.CNV == 'Normal':
            result = '/'.join(sorted(alleles))
        elif r.CNV == 'PartialDeletion':
            if alleles[0] == alleles[1]:
                result = '/'.join(sorted([alleles[0], '*S1']))
            else:
                result = 'Unassigned'
        else:
            result = 'Unassigned'
        return result

    def __init__(self, df):
        self.results = df.apply(self.one_row, axis=1)

class UGT2B17Genotyper:
    """
    Genotyper for UGT2B17.
    """

    def one_row(self, r):
        if r.CNV == 'DeletionHet':
            result = '*1/*2'
        elif r.CNV == 'DeletionHom':
            result = '*2/*2'
        elif r.CNV == 'Normal':
            result = '*1/*1'
        else:
            result = 'Unassigned'
        return result

    def __init__(self, df):
        self.results = df.apply(self.one_row, axis=1)

def call_genotypes(alleles=None, cnv_calls=None):
    """
    Call genotypes for target gene.

    Parameters
    ----------
    alleles : str or pypgx.Archive, optional
        Archive file or object with the semantic type SampleTable[Alleles].
    cnv_calls : str or pypgx.Archive, optional
        Archive file or object with the semantic type SampleTable[CNVCalls].

    Returns
    -------
    pypgx.Archive
        Archive object with the semantic type SampleTable[Genotypes].
    """
    sv_genotypers = {
        'CYP2A6': CYP2A6Genotyper,
        'CYP2B6': CYP2B6Genotyper,
        'CYP2D6': CYP2D6Genotyper,
        'CYP2E1': CYP2E1Genotyper,
        'GSTM1': GSTM1Genotyper,
        'GSTT1': GSTT1Genotyper,
        'SLC22A2': SLC22A2Genotyper,
        'UGT2B15': UGT2B15Genotyper,
        'UGT2B17': UGT2B17Genotyper,
    }

    if isinstance(alleles, str):
        alleles = sdk.Archive.from_file(alleles)

    if alleles is not None:
        alleles.check('SampleTable[Alleles]')

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

    def one_row(r):
        if 'Haplotype1' in r.index:
            r.Haplotype1 = r.Haplotype1.strip(';').split(';')
            r.Haplotype2 = r.Haplotype2.strip(';').split(';')
            d = {}
            for allele in r.VariantData.strip(';').split(';'):
                fields = allele.split(':')
                if 'default' in allele:
                    d[fields[0]] = []
                else:
                    d[fields[0]] = [fields[1].split(','), [float(x) for x in fields[2].split(',')]]
            r.VariantData = d
        return r

    df = df.apply(one_row, axis=1)

    if gene in sv_genotypers:
        df = sv_genotypers[gene](df).results.to_frame()
    else:
        df = SimpleGenotyper(df).results.to_frame()

    df.columns = ['Genotype']

    metadata = {}
    metadata['Gene'] = gene
    metadata['Assembly'] = assembly
    metadata['SemanticType'] = 'SampleTable[Genotypes]'

    return sdk.Archive(metadata, df)
