import os
import time
import pandas as pd
from os.path import join, dirname, abspath

DIR = dirname(abspath(__file__))
GLIPH2_PATH = join(DIR, 'lib')


def GLIPH2(data, outfile=None):
    os.chdir(GLIPH2_PATH)

    data.to_csv('metarepertoire.txt', index=False, header=False, sep='\t')

    print('Clustering {} sequences with GLIPH2.'.format(len(data)))

    # Perform gliph2 algorithm on test sequences
    t0 = time.time()
    os.system('./irtools.centos -c parameters_metarepertoire')
    t1 = time.time()
    t = t1 - t0

    print('Elapsed time: {} seconds.'.format(t))

    # Reformat gliph2 clustering results
    clusters = pd.DataFrame()
    # nodelist = {'CDR3': [], 'cluster': []}
    with open('metarepertoire_output_cluster.txt', 'r') as f:
        results = f.read().splitlines()
    c = 0
    for line in results:
        columns = line.split(' ')
        motif = columns[3]
        cluster = columns[4:]
        if len(cluster) >= 2:
            nodes = pd.DataFrame({'CDR3':cluster})
            nodes['cluster'] = c
            nodes['motif'] = motif
            clusters = clusters.append(nodes)
            c += 1
            
    # for cluster in clusters:
    #     for seq in clusters[cluster]:
    #         nodelist['CDR3'].append(seq)
    #         nodelist['cluster'].append(cluster)
    # nodelist = pd.DataFrame(nodelist)

    if outfile:
        print('Saving output to: \n --> {}'.format(outfile))
        nodelist.to_csv(outfile, sep='\t', index=False)

    return clusters, t
