import os
import random
import pandas as pd

from clustcr.input.vdjdb import parse_vdjdb
from clustcr.input.immuneaccess import parse_immuneaccess
from clustcr.input.tcrex import parse_tcrex
from clustcr.input.airr import parse_airr
from clustcr.input.tenx import parse_10x
from os.path import join, dirname, abspath

DIR = dirname(abspath(__file__))
vdjdb_location = join(DIR, 'vdjdb/vdjdb_full.txt')


def test_cdr3():
    """
    Small data set consisting of 2851 unique CDR3 sequences, curated from a
    subset of the VDJdb.
    This data can be used for testing and benchmarking.
    """
    return vdjdb_beta(q=1, epitopes=False)


def test_epitopes():
    """
    Epitope data corresponding to the sequences in test_cdr3().
    This data can be used for testing and benchmarking.
    """
    return vdjdb_beta(q=1, epitopes=True)


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


def vdjdb_alpha(q=0, epitopes=False):
    vdjdb = parse_vdjdb(vdjdb_location, q=q)
    alpha = vdjdb[['cdr3.alpha', 'v.beta', 'antigen.epitope']].dropna().drop_duplicates()
    alpha = alpha.rename(columns={'cdr3.alpha':'junction_aa',
                                  'v.beta':'v_call',
                                  'antigen.epitope':'epitope'})
    if epitopes:
        return alpha
    else:
        return alpha["junction_aa"].drop_duplicates()

def vdjdb_beta(q=0, epitopes=False):
    vdjdb = parse_vdjdb(vdjdb_location, q=q)
    beta = vdjdb[['cdr3.beta', 'v.beta', 'antigen.epitope']].dropna().drop_duplicates()
    beta = beta.rename(columns={'cdr3.beta':'junction_aa',
                                'v.beta':'v_call',
                                'antigen.epitope':'epitope'})
    if epitopes:
        return beta.reset_index(drop=True)
    else:
        return beta[["junction_aa", "v_call"]].drop_duplicates().reset_index(drop=True)
    
def vdjdb(q=0):
    df = parse_vdjdb(vdjdb_location,q=q)
    df = df.rename(columns={'cdr3.beta':'junction_aa',
                            'v.beta':'v_call',
                            'antigen.epitope':'epitope'})
    return df[["junction_aa", "v_call", "epitope"]].dropna().drop_duplicates()
    
    
def vdjdb_paired(q=0, epitopes=False):
    vdjdb = parse_vdjdb(vdjdb_location, q=q)
    paired = vdjdb[['cdr3.alpha', 'cdr3.beta', 'antigen.epitope']].dropna().drop_duplicates()
    paired.rename(columns={'cdr3.alpha':'CDR3_alpha',
                           'cdr3.beta':'CDR3_beta',
                           'antigen.epitope':'Epitope'},
                  inplace=True)
    if epitopes:
        return paired
    else:
        return paired[['CDR3_alpha', 'CDR3_beta']].drop_duplicates()
    
def vdjdb_gliph2(filename, q=0):
    prepared_data = parse_vdjdb(filename, q=q)
    prepared_data = prepared_data[['cdr3.beta', 'v.beta']].dropna().drop_duplicates()
    prepared_data.rename(columns={'cdr3.beta':'CDR3',
                                  'v.beta':'V'})
    return prepared_data

def vdjdb_tcrdist(filename, q=0):
    prepared_data = parse_vdjdb(filename, q=q)
    prepared_data = prepared_data[['cdr3.beta', 'v.beta']].dropna().drop_duplicates()
    prepared_data.rename(columns={'cdr3.beta':'cdr3_b_aa',
                                  'v.beta':'v_b_gene'})
    return prepared_data
