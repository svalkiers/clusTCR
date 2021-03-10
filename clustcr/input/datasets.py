import os
import random
import pandas as pd

from clustcr.input.vdjdb import vdjdb_to_cdr3list, vdjdb_to_gliph2, vdjdb_to_tcrdist, vdjdb_to_epitopedata
from clustcr.input.immuneaccess import parse_immuneaccess
from clustcr.input.tcrex import parse_tcrex
from clustcr.input.airr import parse_airr
from clustcr.input.tenx import parse_10x
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
    if data_format.lower()=='immuneaccess':
        return parse_immuneaccess(file)
    elif data_format.lower()=='airr':
        return parse_airr(file)
    elif data_format.lower()=='tcrex':
        return parse_tcrex(file)
    elif data_format.lower()=='10x':
        return parse_10x(file)
    else:
        print(f'Unrecognised format: {data_format}')
        

def metarepertoire(directory, data_format, out_format='CDR3', n_sequences=10**6):
    files = os.listdir(directory)
    random.shuffle(files)
    
    if data_format.lower()=='immuneaccess':
        meta = parse_immuneaccess(join(directory, files[0]), out_format=out_format)
    elif data_format.lower()=='airr':
        meta = parse_airr(join(directory, files[0]), out_format=out_format)
    elif data_format.lower()=='tcrex':
        meta = parse_tcrex(join(directory, files[0]), out_format=out_format)

    for file in files[1:]:
        file = join(directory, file)
        if data_format.lower()=='immuneaccess':
            meta = pd.concat([meta, parse_immuneaccess(file, out_format=out_format)])
        elif data_format.lower()=='airr':
            meta = pd.concat([meta, parse_airr(file)])
        elif data_format.lower()=='tcrex':
            meta = pd.concat([meta, parse_tcrex(file)])
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