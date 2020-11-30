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

bm_file = "../data/vdjdb_trb.tsv" # input
bm_output = "../results/bm_metrics.tsv" # output

def read_vdjdb(vdjdb_file, q_score = 0):
    '''
    Reads vdjdb file, keep CDR3 and epitope columns, return entries
    that satisfy the predefined q-score and corresponding CDR3 sequences
    '''        
    vdjdb = pd.read_csv(vdjdb_file, sep="\t")
    vdjdb = vdjdb[vdjdb["Score"]>=q_score]
    
    return vdjdb[["CDR3", "Epitope"]]

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
    bm = read_vdjdb(bm_file)
    cdr3 = set(bm["CDR3"])
    
    # Perform clustering procedure using networkTCR
    Clust = clustering(cdr3)
    edges = Clust.createNetwork()
    nodes = Clust.network_clustering(edges, mcl_hyper=pair)
    
    # Calculate benchmarking metrics
    Benchmark = metrics(nodes, bm)
    retention = Benchmark.retention()
    purity = Benchmark.purity()
    results.append([pair[0], pair[1], retention, purity["Regular"], purity["Baseline"]])
    
    c += 1
    print(c)

# Write output to file
colnames = ["Inflation", "Expansion", "Retention", "Purity_r", "Purity_b"]
output = pd.DataFrame(np.array(results), columns=colnames)
output.to_csv(bm_output, index=False, sep="\t")