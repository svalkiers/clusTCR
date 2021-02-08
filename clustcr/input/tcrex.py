import pandas as pd

def TCRex_to_cdr3list(filename, separator='\t'):
    return pd.read_csv(filename, sep=separator).CDR3_beta.unique()