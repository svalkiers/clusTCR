# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 13:27:39 2020

@author: Sebastiaan Valkiers

Contact: sebastiaan.valkiers@uantwerpen.be
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Results
res = {"ssc":[], "log":[], "size":[]}

# Calculate average information content per amino acid position in cluster
for clust in nodes["cluster"].unique():
    
    sequences = nodes[nodes["cluster"]==clust]["CDR3"].tolist() # sequences of one cluster
    n = len(sequences) # size of cluster
    l = len(sequences[0][1:-1]) # CDR3 length (ignoring pos 1 and -1)
    ic_ssc = [] # small-sample correction
    ic_log = [] # log2(n) correction
    correction = (1/np.log(2))*((20-1)/(2*n)) # small-sample correction
    
    # Make sure to proceed only if all sequences in the cluster have equal length
    if all(len(seq) == len(sequences[0]) for seq in sequences) is False:
        
        # On the rare occasion that a cluster contains sequences of inequal length.
        # Typically, there is/are only one (or very few) sequence(s) that differ from the avg. CDR3 length in the cluster.
        # Therefore, we use the length of the highest proportion of sequences as the standard, and delete all others.
        s = []
        for i in sequences:
            s.append(len(i))
        k = pd.Series(s).value_counts().index[0] # Standard cluster length
        for j in sequences:
            if len(j) != k:
                del sequences[sequences.index(j)] # Delete all sequences that differ from k in length.
    
    else:
        
        for pos in range(1,l+1): # iterate over position 1 to -1 (leaving out first C and last F)
            ic_aa_pos = [] # initiate variable to carry information content per position
            psc = [seq[pos] for seq in sequences] # get first AA of each sequence at position "pos
            
            for aa in pd.Series(psc).value_counts():
                e = (aa/n)*(np.log2(aa/n)) # individual terms of Shannon entropy
                ic_aa_pos.append(e)
            information_content = -sum(ic_aa_pos)
            ic_ssc.append(np.log2(20) - (information_content + correction))
            ic_log.append(information_content / np.log2(n))
        
    res["ssc"].append(np.average(ic_ssc))
    res["log"].append(np.average(ic_log))
    res["size"].append(n)
    
data = pd.DataFrame(res)

fig, ax = plt.subplots(figsize=(12,8))
ax.scatter(data["size"],data["log"])
ax.set_xlabel("cluster size", fontsize=16)
ax.set_ylabel(r"$R_{i}$ per AA position", fontsize=16)
ax.set_title("Cluster variation - VDJdb", fontsize=20)
ax.legend(fontsize=16, loc='lower right')

# for pos in range(1,len(sequences)-1):

    # psc = [seq[1] for seq in tc]
    
    

        # except IndexError:
        #     s = []
        #     for i in X:
        #         s.append(len(i))
        #     k = pd.Series(s).value_counts().index[1]
        #     for j in X:
        #         if len(j) == k:
        #             del X[X.index(j)]
        #     psc = [seq[pos] for seq in X]

# ic = []
# for pos in range(1,len(X[0])-1):
#     psc = [seq[pos] for seq in X]
#     for idx in pd.Series(psc).value_counts():
#         e = (idx/len(psc))*np.log(idx/len(psc))
#         ic.append(e)
#         print(e)