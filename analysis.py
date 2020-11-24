# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 12:17:38 2020

@author: Sebastiaan Valkiers

Contact: sebastiaan.valkiers@uantwerpen.be
"""

from networkTCR import *

data = set(pd.read_csv(r"C:\Users\sebas\Desktop\PhD\Data\Clustering\Benchmarking\nodelist.tsv", sep="\t")["TCR"])

TCR = networkTCR(data)

edges = TCR.createNetwork()
nodes = TCR.network_clustering(edges)

# profile = TCR.calc_profile(nodes)

# consensus_naive = ''.join(profile.idxmax().values)

for cid in nodes["cluster"].unique():
    cluster = nodes[nodes["cluster"]==cid]["TCR"].tolist()
    profile = TCR.calc_profile(cluster)
    print(TCR.clustermotif(profile))