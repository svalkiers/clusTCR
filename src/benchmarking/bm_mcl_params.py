# -*- coding: utf-8 -*-
"""
author: Sebastiaan Valkiers
"""

import pandas as pd
import numpy as np
import time

from networkTCR import Clustering, Metrics

bm_file = "../data/vdjdb_trb.tsv" # input
bm_output = "../results/bm_metrics_t_estimate.tsv" # output

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

# Read input file
bm = read_bm_file(bm_file)
cdr3 = set(bm["CDR3"])

# Generate list of hyperparameter pairs to test
params = []
param_1 = [1.2] 
param_2 = np.round(np.arange(2,10,1),0)
for i in param_1:
    for j in param_2:
        params.append([i, j])

# Scan hyperparameter space
results = []
c = 0
for pair in params:
    
    t0 = time.time()
    # Perform clustering procedure using networkTCR
    print("Clustering...")
    Clust = Clustering(cdr3)
    edges = Clust.createNetwork()
    nodes = Clust.network_clustering(edges, mcl_hyper=pair)
    t1 = time.time()
    
    # Calculate benchmarking metrics
    print("Calculating clustering metrics...")
    Benchmark = Metrics(nodes, bm)
    retention = Benchmark.retention()
    purity = Benchmark.purity()
    cm = Benchmark.calc_confmat()
    consistency = Benchmark.consistency(cm)
    t = t1 - t0
    results.append([t, pair[0], pair[1], retention, purity["True"], purity["Baseline"], consistency["True"], consistency["Baseline"]])
    
    c += 1
    print(c)

# Write output to file
colnames = ["time", "Inflation", "Expansion", "Retention", "Purity_t", "Purity_b", "Consistency_t", "Consistency_b"]
output = pd.DataFrame(np.array(results), columns=colnames)
output.to_csv(bm_output, index=False, sep="\t")