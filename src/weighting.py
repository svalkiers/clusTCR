# -*- coding: utf-8 -*-
"""
author: Sebastiaan Valkiers
"""

import pandas as pd

from networkTCR import Clustering, Dataloader, Metrics

Data = Dataloader('../data/vdjdb_trb.tsv')
vdj = Data.read_bm_file(q = 1)

with open('../data/vdjdb_cdr3_q1.txt', 'r') as f:
    cdr3 = f.read().splitlines()
    indices = [i for i in range(len(cdr3))]
    seqids = pd.DataFrame({"CDR3":cdr3,
                           "ids":indices})

res = {"w":[],
       "d":[],
       "spc":[],
       "retention":[],
       "purity_actual":[],
       "purity_baseline":[],
       "consistency_actual":[],
       "consistency_baseline":[]}

for w in range(1, 6):
    for x in range(1, 11):
        Clust = Clustering(set())
        edges = Clust.import_weighted_network('../data/faissmethod_data/distances.json', n = x, weight = -w)
        nodes = Clust.network_clustering(edgelist = edges)
        nodes.rename(columns={"CDR3":"ids"}, inplace = True)
        results = pd.merge(nodes, seqids, on = 'ids')
        results.drop(columns = ["ids"], inplace = True)
        avg_clust_size = results.shape[0] / results["cluster"].max()
        
        Metr = Metrics(results, vdj)
        ret = Metr.retention()
        pur = Metr.purity()
        con = Metr.consistency()
        res["w"].append(w)
        res["d"].append(x)
        res["spc"].append(avg_clust_size)
        res["retention"].append(ret)
        res["purity_actual"].append(pur['Actual'])
        res["purity_baseline"].append(pur['Baseline'])
        res["consistency_actual"].append(con['Actual'])
        res["consistency_baseline"].append(con['Baseline'])
        
res = pd.DataFrame(res)