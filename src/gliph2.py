import os
import time
import pandas as pd

from load_files.import_functions import path_in_data

def GLIPH2(data):
    
    # Change working directory to gliph2 folder
    main_dir = os.getcwd()
    if main_dir != os.path.abspath('./modules/gliph2/'):
        os.chdir(os.path.abspath('./modules/gliph2/'))
    
    data.to_csv('metarepertoire.txt', index = False, header = False, sep = '\t')
    
    print('Clustering {} sequences with GLIPH2.'.format(len(data)))
    
    # Perform gliph2 algorithm on test sequences
    t0 = time.time()
    os.system('./irtools.centos -c parameters_metarepertoire')
    t1 = time.time()
    t = t1 - t0
    
    print('Elapsed time: {} seconds.'.format(t))
    
    # Reformat gliph2 clustering results
    clusters = {}
    nodelist = {'CDR3':[], 'cluster':[]}
    with open('metarepertoire_output_cluster.txt', 'r') as f:
        results = f.read().splitlines()
    c = 0
    for line in results:
        cluster = line.split(' ')[4:]
        if len(cluster) >= 2:
            clusters[c] = cluster
            c += 1      
    for cluster in clusters:
        for seq in clusters[cluster]:
            nodelist['CDR3'].append(seq)
            nodelist['cluster'].append(cluster)        
    nodelist = pd.DataFrame(nodelist)
    
    # Save output to correct destination
    os.chdir(main_dir)
    outfile = path_in_data('methods/GLIPH2_{}.tsv'.format(len(data)))
    print('Saving output to: \n --> {}'.format(outfile))
    nodelist.to_csv(outfile, sep = '\t', index = False)
    
    return t