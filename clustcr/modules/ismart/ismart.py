import os
import time
import pandas as pd
from os.path import join, dirname, abspath


DIR = dirname(abspath(__file__))
ISMART_PATH = join(DIR, 'lib')


def iSMART(data, outfile=None):
    os.chdir(ISMART_PATH)
    data.to_csv('input.txt', index=False, header=False, sep='\t')

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
    if outfile:
        print('Saving output to: \n --> {}'.format(outfile))
        clusters.to_csv(outfile, sep='\t', index=False)

    return clusters, t
