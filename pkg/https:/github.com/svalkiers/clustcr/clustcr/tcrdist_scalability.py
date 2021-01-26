import time

from load_files.datasets import metarepertoire_tcrdist
from tcrdist.repertoire import TCRrep
from tcrdist.rep_funcs import  compute_pw_sparse_out_of_memory

for i in [50000, 100000, 200000]:
    print('Calculating pw distances of %s sequences.' % (i))
    
    # Sample metarepertoire
    data = metarepertoire_tcrdist(n_sequences=i)
    
    t0 = time.time()
    
    # Start pw distance calculation using tcrdist
    tr = TCRrep(cell_df = data, 
                organism = 'human', 
                chains = ['beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv',
                compute_distances = False,
                store_all_cdr = False)
    
    
    S, fragments = compute_pw_sparse_out_of_memory( tr = tr,
            row_size      = 500,
            pm_processes  = 12,
            pm_pbar       = True,
            max_distance  = 37,
            matrix_name   = 'rw_beta',
            reassemble    = True,
            cleanup       = True)

    t1 = time.time()
    print("Time elapsed: {} seconds".format(t1 - t0))