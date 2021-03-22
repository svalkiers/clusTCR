import pandas as pd

def parser(filename, v_col, cdr3_col, j_col, out_format='CDR3', sep='\t'):
    
    # Read data
    data = pd.read_csv(filename, sep=sep, low_memory=False)[
        [v_col, cdr3_col, j_col]]
    # Rename to ClusTCR format
    data.rename(columns={cdr3_col:"CDR3", v_col:"V", j_col:"J"}, inplace=True)
    # Remove any duplicates and rows with special characters or missing values
    data = data[~data.CDR3.str.contains(r'[^0-9a-zA-Z]')]
    data.drop_duplicates(inplace=True)
    data.dropna(inplace=True)
    if out_format.upper() == 'CDR3':
        return data.CDR3.unique()
    elif out_format.upper() == 'GLIPH2':
        return data
    elif out_format.upper() == 'TCRDIST':
        return data.rename(columns={'CDR3': 'cdr3_b_aa', 'V': 'v_b_gene'})