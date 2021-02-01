import pandas as pd
import json

from .tools import path_in_data


def parse_vdjdb_file(filename, q=0):
    """
    Parse files in the VDJdb format.
    q-score defines the quality of the database entry (3 > 2 > 1 > 0).
    """
    vdjdb = pd.read_csv(path_in_data(filename), sep='\t')
    vdjdb = vdjdb[['CDR3', 'V', 'Score', 'Meta', 'Epitope']]  # Keep necessary columns
    vdjdb = vdjdb[vdjdb['Score'] >= q]  # Quality score cut-off
    vdjdb = vdjdb[~vdjdb.V.str.contains("/")]  # Remove ambiguous entries
    ids = []
    for row in range(len(vdjdb)):
        metadata = json.loads(vdjdb['Meta'].iloc[row])
        ids.append(metadata['subject.id'])
    vdjdb.drop(columns=['Meta', 'Score'], inplace=True)
    vdjdb['subject'] = ids
    vdjdb = vdjdb[vdjdb['subject'].astype(bool)]
    vdjdb = vdjdb[~vdjdb.subject.str.contains('mouse')]
    vdjdb['count'] = [1] * len(vdjdb)
    vdjdb.drop_duplicates(inplace=True)

    return vdjdb


def vdjdb_to_gliph2(filename, q=0):
    prepared_data = parse_vdjdb_file(path_in_data(filename), q=q)
    return prepared_data.drop(columns=['Epitope'], axis=1)


def vdj_to_tcrdist(filename, q=0):
    prepared_data = parse_vdjdb_file(path_in_data(filename), q=q)
    prepared_data.rename(columns={'CDR3': 'cdr3_b_aa', 'V': 'v_b_gene'}, inplace=True)
    return prepared_data.drop(columns=['Epitope'], axis=1)


def vdjdb_to_cdr3list(filename, q=0):
    prepared_data = parse_vdjdb_file(path_in_data(filename), q=q)
    return prepared_data.CDR3


def vdjdb_to_epitopedata(filename, q=0):
    prepared_data = parse_vdjdb_file(path_in_data(filename), q=q)
    return prepared_data.drop(columns=['subject', 'count'], axis=1)
