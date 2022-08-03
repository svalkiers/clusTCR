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
        df = df[df['v_gene'] != 'unresolved']
        df = df.rename(columns={'amino_acid': 'junction_aa',
                                'v_gene': 'v_call'})

    if sample_type == 'v2':
        df = df[df['sequenceStatus'] == 'In']
        df = df[['aminoAcid', 'vGeneName', 'frequencyCount (%)']]
        df = df[df['vGeneName'] != 'unresolved']
        df.rename(columns={'aminoAcid': 'junction_aa',
                           'vGeneName': 'v_call'})

    # Convert Adaptive to IMGT nomenclature
    df['v_call'] = df['v_call'].apply(lambda x: adaptive_to_imgt_human.get(x))
    v_db = imgt_v_genes()
    df = df[df['v_call'].isin(v_db)]

    df.drop_duplicates(inplace=True)
    df['subject'] = [filename.split('/')[-1].replace('.tsv', '')] * len(df)
    df['count'] = df['count'] / 100
    
    if out_format.upper() == 'CDR3':
        return pd.Series(df.junction_aa.unique())
    elif out_format.upper() == 'GLIPH2':
        return df
    elif out_format.upper() == 'TCRDIST':
        return df.rename(columns={'junction_aa': 'cdr3_b_aa', 'v_call': 'v_b_gene'})
