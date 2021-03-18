import pandas as pd

def parse_10X(filename, out_format='CDR3', separator='\t'):
    data = pd.read_csv(filename, sep=separator)
    if out_format=='CDR3':
        return pd.Series(data.Bcdr3.unique())
