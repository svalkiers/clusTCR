from .input import process_csv, path_in_data


def covid19_repertoire():
    filename = path_in_data('covid19_repertoires.tsv')
    return process_csv(filename, cdr3_column_name='CDR3_beta', sep='\t')


def small_vdj():
    filename = path_in_data('VDJ_known_dataset.txt')
    return process_csv(filename, cdr3_column_name='CDR3', epitope_column_name='Epitope', sep='\t')


def vdj():
    filename = path_in_data('VDJdb_human_trb_2020_10_23.tsv')
    return process_csv(filename, cdr3_column_name='CDR3', epitope_column_name='Epitope', sep='\t')


def immune_race():
    files = [
        path_in_data('immune_race_release2.1/peptide-detail-ci.csv'),
        path_in_data('immune_race_release2.1/peptide-detail-cii.csv')
    ]
    return process_csv(files,
                       cdr3_column_name='TCR BioIdentity',
                       cdr3_mapping_func=lambda x: x.split('+')[0],
                       epitope_column_name='Amino Acids')
