import pandas as pd

def parse_airr(filename, out_format='CDR3', separator='\t'):
    data = pd.read_csv(filename, sep=separator)
    # data = data[data["productive"]==True]
    # data = data[data["vj_in_frame"]==True]
    data = data[["aaSeqCDR3", "allVHitsWithScore", "allJHitsWithScore"]]
    data.rename(columns={"aaSeqCDR3":"CDR3", "allVHitsWithScore":"V", "allJHitsWithScore":"J"}, inplace=True)
    data = data[~data.CDR3.str.contains("_")]
    data.drop_duplicates(inplace=True)
    data.dropna(inplace=True)
    if out_format.upper() == 'CDR3':
        return data.CDR3.unique()
    elif out_format.upper() == 'GLIPH2':
        return data
    elif out_format.upper() == 'TCRDIST':
        return data.rename(columns={'CDR3': 'cdr3_b_aa', 'V': 'v_b_gene'})