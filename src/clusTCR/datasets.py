from .import_vdjdb import vdjdb_to_cdr3list, vdjdb_to_gliph2, vdj_to_tcrdist, vdjdb_to_epitopedata
from .import_immuneaccess import construct_metarepertoire, immuneACCESS_to_cdr3list

vdjdb_location = 'vdjdb/vdjdb_trb.tsv'


def test_cdr3():
    """
    Small data set consisting of 2851 unique CDR3 sequences, curated from a
    subset of the VDJdb.
    This data can be used for testing and benchmarking.
    """
    return vdj_cdr3_small()


def test_epitope():
    """
    Epitope data corresponding to the sequences in test_cdr3().
    This data can be used for testing and benchmarking.
    """
    return vdj_epitopes_small()


def vdj_cdr3():
    return vdjdb_to_cdr3list(vdjdb_location)

def vdj_gliph2():
    return vdjdb_to_gliph2(vdjdb_location)

def vdj_epitopes():
    return vdjdb_to_epitopedata(vdjdb_location)

def vdj_cdr3_small(q = 1):
    return vdjdb_to_cdr3list(vdjdb_location, q = q).drop_duplicates()

def vdj_gliph2_small(q = 1):
    return vdjdb_to_gliph2(vdjdb_location, q = q).drop_duplicates()

def vdj_tcrdist_small(q = 1):
    return vdj_to_tcrdist(vdjdb_location, q = q).drop_duplicates()

def vdj_epitopes_small(q = 1):
    return vdjdb_to_epitopedata(vdjdb_location, q = q).drop_duplicates()

def immuneACCESS_cdr3(file):
    return immuneACCESS_to_cdr3list(file)

def metarepertoire_gliph2(directory, n_sequences = 10**6):
    return construct_metarepertoire(directory, n_sequences = n_sequences)

def metarepertoire_tcrdist(directory, n_sequences = 10**6):
    metarepertoire = construct_metarepertoire(directory, n_sequences = n_sequences)
    return metarepertoire.rename(columns = {'CDR3':'cdr3_b_aa', 'V':'v_b_gene'})

def metarepertoire_cdr3(directory, n_sequences = 10**6):
    return construct_metarepertoire(directory, n_sequences = n_sequences).CDR3