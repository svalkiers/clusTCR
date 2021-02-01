import pandas as pd
import os
import random
from os.path import join

from .tools import path_in_data, imgt_v_genes
from clustcr.input.adaptive_to_imgt import adaptive_to_imgt_human


def parse_immuneACCESS(filename, separator='\t'):
    """
    Parse data in the immuneACCESS format.
    """
    df = pd.read_csv(filename, sep=separator)
    if 'amino_acid' in df.columns:
        sample_type = 'v1'
    elif 'aminoAcid' in df.columns:
        sample_type = 'v2'
    else:
        raise RuntimeError('immuneACCESS invalid input format')

    if sample_type == 'v1':
        df = pd.read_csv(filename, sep=separator)
        df = df[df['frame_type'] == 'In']
        df = df[['amino_acid', 'v_gene', 'productive_frequency']]
        df = df[df['vGeneName'] != 'unresolved']
        df.rename(columns={'amino_acid': 'CDR3',
                           'v_gene': 'V',
                           'productive_frequency': 'count'},
                  inplace=True)

    if sample_type == 'v2':
        df = df[df['sequenceStatus'] == 'In']
        df = df[['aminoAcid', 'vGeneName', 'frequencyCount (%)']]
        df = df[df['vGeneName'] != 'unresolved']
        df.rename(columns={'aminoAcid': 'CDR3',
                           'vGeneName': 'V',
                           'frequencyCount (%)': 'count'},
                  inplace=True)

    # Convert Adaptive to IMGT nomenclature
    df['V'] = df['V'].apply(lambda x: adaptive_to_imgt_human.get(x))
    v_db = imgt_v_genes()
    df = df[df['V'].isin(v_db)]

    df.drop_duplicates(inplace=True)
    df['subject'] = [filename.split('/')[-1].replace('.tsv', '')] * len(df)
    df['count'] = df['count'] / 100

    return df


def immuneACCESS_to_gliph2(filename):
    return parse_immuneACCESS(path_in_data(filename))


def immuneACCESS_to_cdr3list(filename):
    prepared_data = parse_immuneACCESS(path_in_data(filename))
    return prepared_data.CDR3


def construct_metarepertoire(directory, n_sequences=10 ** 6):
    files = os.listdir(directory)
    random.shuffle(files)
    meta = parse_immuneACCESS(join(directory, files[0]))

    for file in files[1:]:
        file = join(directory, file)
        meta = pd.concat([meta, parse_immuneACCESS(file)])
        meta.drop_duplicates(inplace=True)
        if len(meta) > n_sequences:
            return meta.sample(n_sequences)

    print(f'Metarepertoire: less sequences found than wanted ({len(meta)} vs {n_sequences})')
    return meta
