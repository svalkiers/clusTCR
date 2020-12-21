# -*- coding: utf-8 -*-
"""
author: Sebastiaan Valkiers
"""

import numpy as np

from networkTCR import Dataloader, Clustering, Features, Metrics
from clusterAnalyzer import ClusterAnalysis

# Load data
Data = Dataloader("../data/vdjdb_trb.tsv")
vdjdb = Data.read_bm_file(q=0)

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
cm_t, cm_b = Metr.calc_confmat()
labels_t = []
labels_b = []
for i in cm_t:
    labels_t.append(cm_t[i].max() / np.sum(cm_t[i]))
    labels_b.append(cm_b[i].max() / np.sum(cm_b[i]))

# Predict purity
Analyzer = ClusterAnalysis(features, labels_t, labels_b)
Analyzer.plot_purity_distribution()
Analyzer.discretize_labels()
Analyzer.prep_data()
Analyzer.cluster_quality_classifier(feat_importances = True)
Analyzer.perform_PCA()