import pandas as pd

def parse_airr(filename, separator='\t'):
    data = pd.read_csv(filename, sep=separator)
    data = data[data["productive"]==True]
    data = data[data["vj_in_frame"]==True]
    data = data[["sequence_aa", "v_call", "j_call"]]
    data.rename(columns={"sequence_aa":"CDR3", "v_call":"V", "j_call":"J"}, inplace=True)
    data.drop_duplicates(inplace=True)
    return data

def airr_to_cdr3list(filename):
    return parse_airr(filename).CDR3.unique()