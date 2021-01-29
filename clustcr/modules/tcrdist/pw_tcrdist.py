import time

from tcrdist.repertoire import TCRrep
from tcrdist.rep_funcs import  compute_pw_sparse_out_of_memory


def tcrdist(data):
    
    print('pw dist calculations for {} sequences with tcrdist.'.format(len(data)))
    
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
            pm_processes  = 1,
            pm_pbar       = True,
            max_distance  = 37,
            matrix_name   = 'rw_beta',
            reassemble    = True,
            cleanup       = True)

    t1 = time.time()
    t = t1 - t0
    
    print('Elapsed time: {} seconds.'.format(t1 - t0))
    
    return S, t