import pandas as pd
import os
import random

from .import_functions import path_in_data


def parse_immuneACCESS(filename, separator = '\t'):
    '''
    Parse data in the immuneACCESS format.
    '''
    df = pd.read_csv(filename, sep = separator)
    df = df[df['sequenceStatus'] == 'In']
    df = df[['aminoAcid', 'vGeneName']]
    df = df[df['vGeneName'] != 'unresolved']
    df.rename(columns = {'aminoAcid' : 'CDR3',
                         'vGeneName' : 'V'},
              inplace = True)
    df.drop_duplicates(inplace = True)
    df['subject'] = [filename.split('/')[-1].replace('.tsv', '')] * len(df)
    df['count'] = [1] * len(df)
    
    return df

def immuneACCESS_to_gliph2(filename):
    return parse_immuneACCESS(path_in_data(filename))

def immuneACCESS_to_cdr3list(filename):
    prepared_data = parse_immuneACCESS(path_in_data(filename))
    return prepared_data.CDR3

# def construct_metarepertoire(directory, n_sequences = 10**6):
#     initiate = True
#     folder = path_in_data(directory)
#     for repertoire in os.listdir(folder):
#         if initiate:
#             meta = parse_immuneACCESS(folder + repertoire)
#             meta.drop_duplicates(inplace = True)
#             initiate = False
#             if len(meta) >= n_sequences:
#                 meta = meta.sample(n_sequences)
#                 break
#         else:
#             rep = parse_immuneACCESS(folder + repertoire)
#             meta = pd.concat([meta, rep])
#             meta.drop_duplicates(inplace = True)
#             if len(meta) >= n_sequences:
#                 meta = meta.sample(n_sequences)
#                 break
#     return meta

def construct_metarepertoire(directory, n_sequences = 10**6):

    folder = path_in_data(directory)
    meta = parse_immuneACCESS(folder + random.choice(os.listdir(folder)))
    meta.drop_duplicates(inplace = True)
    
    while len(meta) <= n_sequences:
        rep = parse_immuneACCESS(folder + random.choice(os.listdir(folder)))
        meta = pd.concat([meta, rep])
        meta.drop_duplicates(inplace = True)
        
    meta = meta.sample(n_sequences)
    
    return meta