import time
import pandas as pd
from os.path import join, dirname, abspath, exists
import multiprocessing
import parmap
import random
from os import mkdir
from subprocess import call
from shutil import rmtree


DIR = dirname(abspath(__file__))
ISMART_PATH = join(DIR, 'lib')
ISMART_EXEC = join(ISMART_PATH, 'iSMARTf3.py')


def iSMART(data, outfile=None):
    tmp_directory = join(DIR, 'ismart_tmp' + str(random.randint(0, 10 ** 8)))
    if exists(tmp_directory):
        rmtree(tmp_directory)
    mkdir(tmp_directory)

    input_file = join(tmp_directory, 'input.txt')
    data.to_csv(input_file, index=False, header=False, sep='\t')

    print('Clustering {} sequences with iSMART.'.format(len(data)))

    # Perform gliph2 algorithm on test sequences
    t0 = time.time()
    # os.system('python iSMARTf3.py -f input.txt -v False')
    call(f'python {ISMART_EXEC} -f {input_file} -o {tmp_directory} -v False', shell=True)
    t1 = time.time()
    t = t1 - t0

    print('Elapsed time: {} seconds.'.format(t))

    with open(join(tmp_directory, 'input_clustered_v3.txt'), 'r') as f:
        clusters = f.read().splitlines()[3:]
        clusters = pd.DataFrame([x.split('\t') for x in clusters], columns=['CDR3', 'cluster'])
        clusters.cluster = pd.to_numeric(clusters.cluster, errors='coerce')

    # Save output to correct destination
    if outfile:
        print('Saving output to: \n --> {}'.format(outfile))
        clusters.to_csv(outfile, sep='\t', index=False)

    rmtree(tmp_directory)

    return clusters, t


def iSMART_from_preclusters(preclusters, n_cpus):
    with multiprocessing.Pool(n_cpus) as pool:
        clusters = parmap.map(iSMART,
                              preclusters,
                              pm_parallel=True,
                              pm_pool=pool)
    for c in range(len(clusters)):
        clusters[c] = clusters[c][0]
        if c != 0:
            clusters[c]['cluster'] += clusters[c - 1]['cluster'].max() + 1
    return pd.concat(clusters, ignore_index=True)

