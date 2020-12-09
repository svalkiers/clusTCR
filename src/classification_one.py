# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 15:26:12 2020

@author: Sebastiaan Valkiers

Contact: sebastiaan.valkiers@uantwerpen.be
"""

import numpy as np
import pandas as pd

from networkTCR import Dataloader, Clustering, Features, Metrics

# Load data
Data = Dataloader("../data/vdjdb_trb.tsv")
vdjdb = Data.read_bm_file(q=1)

# Perform clustering
Clust = Clustering(set(vdjdb["CDR3"]))
edges = Clust.create_network()
nodes = Clust.network_clustering(edges)

# Calculate features of clusters
Feat = Features(nodes)
aavar = Feat.calc_variation()
pchem = Feat.calc_physchem()
pgen = Feat.calc_pgen()
features = Feat.combine(aavar, pchem, pgen) # combines results into joint DataFrame

# Generate labels (purity)
Metr = Metrics(nodes, vdjdb)
cm = Metr.calc_confmat()[0]
labels = []
for i in cm:
    labels.append(cm[i].max() / np.sum(cm[i]))