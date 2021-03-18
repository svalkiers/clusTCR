import pandas as pd

def parse_TCRex(filename, out_format='CDR3', separator='\t'):
    data = pd.read_csv(filename, sep=separator)
    data.drop(columns=['TRBJ_gene'], inplace=True)
    if out_format.upper() == 'CDR3':
        return data.CDR3_beta.unique()
    elif out_format.upper() == 'GLIPH2':
        return data.rename(columns={'CDR3_beta': 'CDR3', 'TRBV_gene': 'V'})
    elif out_format.upper() == 'TCRDIST':
        return data.rename(columns={'CDR3': 'cdr3_b_aa', 'V': 'v_b_gene'})
