import os
import time
import pandas as pd

from load_files.import_functions import path_in_data

def iSMART(data):
    
    # Change working directory to gliph2 folder
    main_dir = os.getcwd()
    if main_dir != os.path.abspath('./modules/ismart/'):
        os.chdir(os.path.abspath('./modules/ismart/'))
    
    data.to_csv('input.txt', index = False, header = False, sep = '\t')
    
    print('Clustering {} sequences with iSMART.'.format(len(data)))
    
    # Perform gliph2 algorithm on test sequences
    t0 = time.time()
    os.system('python iSMARTf3.py -f input.txt -v False')
    t1 = time.time()
    t = t1 - t0
    
    print('Elapsed time: {} seconds.'.format(t))
    
    with open('input_clustered_v3.txt', 'r') as f:
        clusters = f.read().splitlines()[3:]
        clusters = pd.DataFrame([x.split('\t') for x in clusters], columns=['CDR3', 'cluster'])
    
    # Save output to correct destination
    os.chdir(main_dir)
    outfile = path_in_data('methods/iSMART_{}.tsv'.format(len(data)))
    print('Saving output to: \n --> {}'.format(outfile))
    clusters.to_csv(outfile, sep = '\t', index = False)
    
    return t
    