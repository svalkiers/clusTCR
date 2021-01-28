import pandas as pd
from os.path import join, dirname, abspath

DIR = dirname(abspath(__file__))

mapping = pd.read_csv(join(DIR, 'adaptive_imgt_mapping.csv'))

adaptive_to_imgt_human = mapping.loc[mapping['species'] == 'human'].set_index('adaptive')['imgt'].fillna('NA').to_dict()
