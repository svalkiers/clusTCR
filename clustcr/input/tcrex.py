import pandas as pd

def parse_tcrex(filename, out_format='CDR3', separator='\t'):
    data = pd.read_csv(filename, sep=separator)
    data.drop(columns=['TRBJ_gene'], inplace=True)
    data = data.rename(columns={"TRBV_gene":"v_call",
                                "CDR3_beta":"junction_aa"})
    if out_format.upper() == 'CDR3':
        return pd.Series(data.junction_aa.unique())
    elif out_format.upper() == 'GLIPH2':
        return data.rename(columns={'junction_aa':'CDR3', 
                                    'v_call':'V'})
    elif out_format.upper() == 'TCRDIST':
        return data.rename(columns={'junction_aa':'cdr3_b_aa',
                                    'v_call':'v_b_gene'})
