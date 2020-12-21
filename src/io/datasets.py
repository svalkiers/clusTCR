from import_functions import path_in_data, from_tcrdata


def vdj_cdr3():
    filename = path_in_data('vdjdb_trb.tsv')
    return from_tcrdata(filename)

def vdj_with_epitopes():
    filename = path_in_data('vdjdb_trb.tsv')
    return from_tcrdata(filename, epitopes = True)

def vdj_small_cdr3():
    filename = path_in_data('vdjdb_trb.tsv')
    return from_tcrdata(filename, q = 1)

def vdj_small_with_epitopes():
    filename = path_in_data('vdjdb_trb.tsv')
    return from_tcrdata(filename, epitopes = True, q = 1)