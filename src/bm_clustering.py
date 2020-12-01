# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 12:17:38 2020

@author: Sebastiaan Valkiers

Contact: sebastiaan.valkiers@uantwerpen.be
"""

import pandas as pd
import numpy as np

from networkTCR import clustering
from networkTCR import metrics

bm_file = "../data/dash.tsv" # input
bm_output = "../results/bm_metrics_dash.tsv" # output

def read_bm_file(file, q_score = None):
    '''
    Reads benchmark file, keep CDR3 and epitope columns, return entries
    that satisfy the predefined q-score and corresponding CDR3 sequences
    if using VDJdb.
    '''        
    bm = pd.read_csv(file, sep="\t")
    if q_score is not None:
        bm = bm[bm["Score"]>=q_score]
    
    return bm[["CDR3", "Epitope"]]

# Generate list of hyperparameter pairs to test
params = []
param_1 = np.round(np.arange(1.1,3.1,0.1),1)
param_2 = np.round(np.arange(2,6,1),0)
for i in param_1:
    for j in param_2:
        params.append([i, j])

# Scan hyperparameter space
results = []
c = 0
for pair in params:

    # Read input file
    bm = read_bm_file(bm_file)
    cdr3 = set(bm["CDR3"])
    
    # Perform clustering procedure using networkTCR
    print("Clustering...")
    Clust = clustering(cdr3)
    edges = Clust.createNetwork()
    nodes = Clust.network_clustering(edges, mcl_hyper=pair)
    
    # Calculate benchmarking metrics
    print("Calculating clustering metrics...")
    Benchmark = metrics(nodes, bm)
    retention = Benchmark.retention()
    purity = Benchmark.purity()
    cm = Benchmark.calc_confmat()
    consistency = Benchmark.consistency(cm)
    results.append([pair[0], pair[1], retention, purity["Regular"], purity["Baseline"], consistency["Regular"], consistency["Baseline"]])
    
    c += 1
    print(c)

# Write output to file
colnames = ["Inflation", "Expansion", "Retention", "Purity_r", "Purity_b", "Consistency_r", "Consistency_b"]
output = pd.DataFrame(np.array(results), columns=colnames)
output.to_csv(bm_output, index=False, sep="\t")