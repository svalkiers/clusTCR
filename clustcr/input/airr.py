import pandas as pd

def parse_airr(filename, out_format='CDR3', separator='\t'):
    data = pd.read_csv(filename, sep=separator)
    data = data[data["productive"]==True]
    data = data[data["vj_in_frame"]==True]
    data = data[["sequence_aa", "v_call", "j_call"]]
    data.rename(columns={"sequence_aa":"CDR3", "v_call":"V", "j_call":"J"}, inplace=True)
    data.drop_duplicates(inplace=True)
    if out_format.upper() == 'CDR3':
        return data.CDR3.unique()
    elif out_format.upper() == 'GLIPH2':
        return data
    elif out_format.upper() == 'TCRDIST':
        return data.rename(columns={'CDR3': 'cdr3_b_aa', 'V': 'v_b_gene'})