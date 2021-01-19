import time
import pandas as pd

from load_files.datasets import metarepertoire_gliph2, metarepertoire_cdr3, metarepertoire_tcrdist
from gliph2 import GLIPH2
from ismart import iSMART
from pw_tcrdist import tcrdist
from clusTCR import Clustering

# Different repertoire sizes for scalability testing
repertoire_sizes = [100000, 200000]

# Keep track of speed at each size
time_records = {'GLIPH2' : [],
                'iSMART' : [],
                'Two-step' : [],
                'tcrdist' : [],
                'n_sequences' : []}

# Perform clustering at different repertoire sizes
# with different SoA clustering methods
for s in repertoire_sizes:
    
    # Get data
    rep_gliph = metarepertoire_gliph2(s)
    rep_td = metarepertoire_tcrdist(s)
    rep_cdr3 = metarepertoire_cdr3(s)
    
    # Two-step clustering
    t0 = time.time()
    print('Clustering {} sequences with Two-step.'.format(len(rep_cdr3)))
    clstr = Clustering(rep_cdr3, method = 'two-step', n_cpus=8)
    output = clstr.cluster_sequences()
    t1 = time.time()
    time_records['Two-step'].append(t1 - t0)
    print('Elapsed time: {} seconds.'.format(t1 - t0))
    
    # GLIPH2 clustering
    try:
        t_gliph = GLIPH2(rep_gliph)
        time_records['GLIPH2'].append(t_gliph)
    except:
        time_records['GLIPH2'].append(None)
    
    # iSMART clustering
    try:
        t_ismart = iSMART(rep_cdr3)
        time_records['iSMART'].append(t_ismart)
    except:
        time_records['iSMART'].append(None)
        
    # tcrdist pw distance calculation (no clustering)
    t_td = tcrdist(rep_td)
    time_records['tcrdist'].append(t_td)
    
    # Prepare results file
    time_records['n_sequences'].append(s)
    speed = pd.DataFrame(time_records)
    speed.to_csv('../results/method_comparison_speed_small_2.tsv', sep = '\t', index = False)