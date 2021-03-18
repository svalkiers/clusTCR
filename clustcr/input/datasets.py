import os
import random
import pandas as pd

from clustcr.input.VDJdb import vdjdb_to_cdr3list, vdjdb_to_gliph2, vdjdb_to_tcrdist, vdjdb_to_epitopedata
from clustcr.input.AIRR import parse_AIRR
from clustcr.input.immuneACCESS import parse_immuneACCESS
from clustcr.input.TCRex import parse_TCRex
from clustcr.input.MiXCR import parse_MiXCR
from clustcr.input.TenX import parse_10X
from os.path import join, dirname, abspath

DIR = dirname(abspath(__file__))
vdjdb_location = join(DIR, 'vdjdb_trb.tsv')


def test_cdr3():
    """
    Small data set consisting of 2851 unique CDR3 sequences, curated from a
    subset of the VDJdb.
    This data can be used for testing and benchmarking.
    """
    return vdjdb_cdr3_small()


def test_epitopes():
    """
    Epitope data corresponding to the sequences in test_cdr3().
    This data can be used for testing and benchmarking.
    """
    return vdjdb_epitopes_small()


def read_cdr3(file, data_format):
    """
    Import function to read and parse rep-seq data of various formats.
    """
    if data_format.lower()=='airr':
        return parse_AIRR(file)
    elif data_format.lower()=='immuneaccess':
        return parse_immuneACCESS(file)
    elif data_format.lower()=='mixcr':
        return parse_MiXCR(file)
    elif data_format.lower()=='tcrex':
        return parse_TCRex(file)
    elif data_format.lower()=='10x':
        return parse_10X(file)
    else:
        print(f'Unrecognised format: {data_format}')
        

def metarepertoire(directory, data_format, files=None, out_format='CDR3', n_sequences=10**6):
    
    if files is None:
        print("Randomly sampling from directory:\n%s" % directory)
        files = os.listdir(directory)
        random.shuffle(files)
    else:
        print("Randomly sampling from provided file list...")
    
    meta = pd.DataFrame()

    for file in files:
        file = join(directory, file)
        if data_format.lower()=='airr':
            meta = pd.concat([meta, parse_AIRR(file, out_format=out_format)])
        elif data_format.lower()=='immuneaccess':
            meta = pd.concat([meta, parse_immuneACCESS(file, out_format=out_format)])
        elif data_format.lower()=='mixcr':
            meta = pd.concat([meta, parse_MiXCR(file)])
        elif data_format.lower()=='tcrex':
            meta = pd.concat([meta, parse_TCRex(file)])
        elif data_format.lower()=='10x':
            meta = pd.concat([meta, parse_10X(file)])
        meta.drop_duplicates(inplace=True)
        if len(meta) > n_sequences:
            return meta.sample(n_sequences)

    print(f'Metarepertoire: less sequences found than desired ({len(meta)} vs {n_sequences})')
    return meta


def vdjdb_cdr3():
    return vdjdb_to_cdr3list(vdjdb_location)


def vdjdb_gliph2():
    return vdjdb_to_gliph2(vdjdb_location)


def vdjdb_epitopes():
    return vdjdb_to_epitopedata(vdjdb_location)


def vdjdb_cdr3_small(q=1):
    return vdjdb_to_cdr3list(vdjdb_location, q=q).drop_duplicates()


def vdjdb_gliph2_small(q=1):
    return vdjdb_to_gliph2(vdjdb_location, q=q).drop_duplicates()


def vdjdb_tcrdist_small(q=1):
    return vdjdb_to_tcrdist(vdjdb_location, q=q).drop_duplicates()


def vdjdb_epitopes_small(q=1):
    return vdjdb_to_epitopedata(vdjdb_location, q=q).drop_duplicates()