import pandas as pd

from .tools import imgt_v_genes
from .adaptive_to_imgt import adaptive_to_imgt_human


def parse_immuneaccess(filename, out_format='CDR3', separator='\t'):
    """
    Parse data in the immuneACCESS format.
    """
    df = pd.read_csv(filename, sep=separator)
    if 'amino_acid' in df.columns:
        sample_type = 'v1'
    elif 'aminoAcid' in df.columns:
        sample_type = 'v2'
    else:
        raise RuntimeError('immuneACCESS invalid input format')

    if sample_type == 'v1':
        df = pd.read_csv(filename, sep=separator)
        df = df[df['frame_type'] == 'In']
        df = df[['amino_acid', 'v_gene', 'productive_frequency']]
        df = df[df['vGeneName'] != 'unresolved']
        df.rename(columns={'amino_acid': 'CDR3',
                           'v_gene': 'V',
                           'productive_frequency': 'count'},
                  inplace=True)

    if sample_type == 'v2':
        df = df[df['sequenceStatus'] == 'In']
        df = df[['aminoAcid', 'vGeneName', 'frequencyCount (%)']]
        df = df[df['vGeneName'] != 'unresolved']
        df.rename(columns={'aminoAcid': 'CDR3',
                           'vGeneName': 'V',
                           'frequencyCount (%)': 'count'},
                  inplace=True)

    # Convert Adaptive to IMGT nomenclature
    df['V'] = df['V'].apply(lambda x: adaptive_to_imgt_human.get(x))
    v_db = imgt_v_genes()
    df = df[df['V'].isin(v_db)]

    df.drop_duplicates(inplace=True)
    df['subject'] = [filename.split('/')[-1].replace('.tsv', '')] * len(df)
    df['count'] = df['count'] / 100
    
    if out_format.upper() == 'CDR3':
        return pd.Series(df.CDR3.unique())
    elif out_format.upper() == 'GLIPH2':
        return df
    elif out_format.upper() == 'TCRDIST':
        return df.rename(columns={'CDR3': 'cdr3_b_aa', 'V': 'v_b_gene'})