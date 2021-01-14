from .import_vdjdb import vdjdb_to_cdr3list, vdjdb_to_gliph2, vdjdb_to_epitopedata
from .import_immuneaccess import construct_metarepertoire, immuneACCESS_to_cdr3list


def vdj_cdr3():
    return vdjdb_to_cdr3list('vdjdb/vdjdb_trb.tsv')

def vdj_gliph2():
    return vdjdb_to_gliph2('vdjdb/vdjdb_trb.tsv')

def vdj_epitopes():
    return vdjdb_to_epitopedata('vdjdb/vdjdb_trb.tsv')

def vdj_cdr3_small(q = 1):
    return vdjdb_to_cdr3list('vdjdb/vdjdb_trb.tsv', q = q)

def vdj_gliph2_small(q = 1):
    return vdjdb_to_gliph2('vdjdb/vdjdb_trb.tsv', q = q)

def vdj_epitopes_small(q = 1):
    return vdjdb_to_epitopedata('vdjdb/vdjdb_trb.tsv', q = q)

def immuneACCESS_cdr3(file):
    return immuneACCESS_to_cdr3list(file)

def metarepertoire_gliph2(n_sequences = 10**6):
    return construct_metarepertoire('emerson_cohort_1/', n_sequences = n_sequences)

def metarepertoire_cdr3(n_sequences = 10**6):
    metarepertoire = construct_metarepertoire('emerson_cohort_1/', n_sequences = n_sequences)
    return metarepertoire.CDR3