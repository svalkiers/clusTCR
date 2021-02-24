import os
import time
import pandas as pd
import multiprocessing
import parmap
from os.path import join, dirname, abspath
from subprocess import call

DIR = dirname(abspath(__file__))
GLIPH2_PATH = join(DIR, 'lib')


def GLIPH2(data, outfile=None):
    os.chdir(GLIPH2_PATH)

    data.to_csv('metarepertoire.txt', index=False, header=False, sep='\t')

    print('Clustering {} sequences with GLIPH2.'.format(len(data)))

    # Perform gliph2 algorithm on test sequences
    t0 = time.time()
    call('./irtools.centos -c parameters_metarepertoire', shell=True)
    t1 = time.time()
    t = t1 - t0

    print('Elapsed time: {} seconds.'.format(t))

    # Reformat gliph2 clustering results
    clusters = {}
    nodelist = {'CDR3': [], 'cluster': []}
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

    if outfile:
        print('Saving output to: \n --> {}'.format(outfile))
        nodelist.to_csv(outfile, sep='\t', index=False)

    return nodelist, t


def GLIPH2_from_preclusters(preclusters, n_cpus):
    with multiprocessing.Pool(n_cpus) as pool:
        clusters = parmap.map(GLIPH2,
                              preclusters,
                              pm_parallel=True,
                              pm_pool=pool)
    for c in range(len(clusters)):
        clusters[c] = clusters[c][0]
        if c != 0:
            clusters[c]['cluster'] += clusters[c - 1]['cluster'].max() + 1
    return pd.concat(clusters, ignore_index=True)


