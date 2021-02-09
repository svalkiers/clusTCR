from clustcr.input.vdjdb import vdjdb_to_cdr3list, vdjdb_to_gliph2, vdjdb_to_tcrdist, vdjdb_to_epitopedata
from clustcr.input.immuneaccess import construct_metarepertoire, immuneACCESS_to_cdr3list
from clustcr.input.tcrex import TCRex_to_cdr3list
from clustcr.input.airr import airr_to_cdr3list
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


def read_cdr3(file, format):
    if format.lower()=='immuneaccess':
        return immuneACCESS_to_cdr3list(file)
    elif format.lower()=='airr':
        return airr_to_cdr3list(file)
    elif format.lower()=='tcrex':
        return TCRex_to_cdr3list(file)
    else:
        print('Unrecognised format: %s' % (format))


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


def metarepertoire_gliph2(location, n_sequences=10 ** 6):
    return construct_metarepertoire(location, n_sequences=n_sequences)


def metarepertoire_tcrdist(location, n_sequences=10 ** 6):
    metarepertoire = construct_metarepertoire(location, n_sequences=n_sequences)
    return metarepertoire.rename(columns={'CDR3': 'cdr3_b_aa', 'V': 'v_b_gene'})


def metarepertoire_cdr3(location, n_sequences=10 ** 6):
    return construct_metarepertoire(location, n_sequences=n_sequences).CDR3
