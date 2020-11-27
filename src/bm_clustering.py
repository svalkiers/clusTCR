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

bm_file = "../data/vdjdb_trb.tsv"
bm_output = "../results/bm_metrics.tsv"

def read_vdjdb(vdjdb_file, q_score = 0):
    '''
    Reads vdjdb file, keep CDR3 and epitope columns, return entries
    that satisfy the predefined q-score and corresponding CDR3 sequences
    '''        
    vdjdb = pd.read_csv(vdjdb_file, sep="\t")
    vdjdb = vdjdb[vdjdb["Score"]>=q_score]
    
    return vdjdb[["CDR3", "Epitope"]]

# Read input file
bm = read_vdjdb(bm_file)
cdr3 = set(bm["CDR3"])

# Perform clustering procedure using networkTCR
Clust = clustering(cdr3)
edges = Clust.createNetwork()
nodes = Clust.network_clustering(edges)

# Calculate benchmarking metrics
Benchmark = metrics(nodes, bm)
retention = Benchmark.retention()
purity = Benchmark.purity()
results = [1.2, 2, retention, purity["Regular"], purity["Baseline"]]

# Write output to file
colnames = ["Inflation", "Expansion", "Retention", "Purity_r", "Purity_b"]
output = pd.DataFrame(np.array(results).reshape(-1,len(results)), columns=colnames)
output.to_csv(bm_output, index=False, sep="\t")