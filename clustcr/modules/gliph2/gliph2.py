import os
import time
import pandas as pd
import multiprocessing
import parmap
from os.path import join, dirname, abspath, exists
from subprocess import call
from shutil import rmtree
from os import mkdir
import random

DIR = dirname(abspath(__file__))
GLIPH2_PATH = join(DIR, 'lib')
GLIPH2_EXEC = join(GLIPH2_PATH, 'irtools.centos')


def GLIPH2(data, outfile=None):
    tmp_directory = join(DIR, 'gliph2_tmp' + str(random.randint(0, 10 ** 8)))
    if exists(tmp_directory):
        rmtree(tmp_directory)

    print('Clustering {} sequences with GLIPH2.'.format(len(data)))
    call(f'cp -r {GLIPH2_PATH} {tmp_directory}', shell=True)
    call(f'chmod a+x {join(tmp_directory, "irtools.centos")}', shell=True)

    input_file = join(tmp_directory, 'metarepertoire.txt')
    data.to_csv(input_file, index=False, header=False, sep='\t')

    with open(join(GLIPH2_PATH, 'parameters_metarepertoire')) as f:
        parameters = f.read()
        parameters = parameters.replace('ref_CD48_v2.0.fa', join(DIR, 'ref/ref_CD48_v2.0.fa'))
        parameters = parameters.replace('ref_V_CD48_v2.0.txt', join(DIR, 'ref/ref_V_CD48_v2.0.txt'))
        parameters = parameters.replace('ref_L_CD48_v2.0.txt', join(DIR, 'ref/ref_L_CD48_v2.0.txt'))

    parameters_file = join(tmp_directory, 'parameters_metarepertoire')
    with open(parameters_file, 'w') as f:
        f.write(parameters)

    # Perform gliph2 algorithm on test sequences
    t0 = time.time()
    call(f'cd {tmp_directory} && ./irtools.centos -c {parameters_file}', shell=True)
    t1 = time.time()
    t = t1 - t0

    print('Elapsed time: {} seconds.'.format(t))

    # Reformat gliph2 clustering results
    clusters = {}
    nodelist = {'CDR3': [], 'cluster': []}
    with open(join(tmp_directory, 'metarepertoire_output_cluster.txt'), 'r') as f:
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

    rmtree(tmp_directory)
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


