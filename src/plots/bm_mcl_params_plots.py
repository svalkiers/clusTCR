# -*- coding: utf-8 -*-
""" 
@author: Sebastiaan Valkiers
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib import cm

# Adjust these variables if other metrics/datasets are available
metrics = ["Purity", "Consistency"]
datasets = ["VDJdb", "Dash"]

for dataset in datasets:
    
    # Prepare data for plotting
    infile = "../results/bm_metrics_" + dataset + ".tsv"
    df = pd.read_csv(infile, sep="\t")
    hyperparameters = df["Expansion"].unique()
    colors = cm.get_cmap('viridis', len(hyperparameters))
    
    for metric in metrics:
        
        # Make new plot for every metric
        fig, ax = plt.subplots(figsize=[12,8])
        
        for hp, col in zip(hyperparameters, colors.colors):
            df_hp = df[df["Expansion"]==int(hp)]
            ax.plot(df_hp["Inflation"], df_hp[metric + "_r"], c=col, label="True. Expansion = " + str(int(hp))) # True
            ax.plot(df_hp["Inflation"], df_hp[metric + "_b"], c=col, label="Base. Expansion = " + str(int(hp)), linestyle="dotted") # Baseline
            
        ax.set_ylabel(metric, fontsize=16)
        ax.set_xlabel("Inflation", fontsize=16)
        ax.set_title("MCL hyperparameters versus " + metric.lower() + " - " + dataset, fontsize=18)
        ax.legend()
        

# Summary plot comparing different q-score cutoffs (VDJdb)
summ = pd.read_csv("../results/summary_VDJdb.tsv", sep="\t")
summ = summ[summ["d"]==1]

fig, axs = plt.subplots(2, 2, figsize=[12,8])

labels = [0,1,2]
x = np.arange(len(labels))
width = 0.35

# Number
axs[0,0].bar(x=x, height=summ["nodes"], tick_label=labels)
axs[0,0].set_title("Number of nodes in network", fontsize=18)
axs[0,0].set_xlabel("VDJdb q-score", fontsize=14)
axs[0,0].set_ylabel("# nodes", fontsize=14)

# Retention
axs[0,1].bar(x=x, height=summ["retention"], tick_label=labels)
axs[0,1].set_title("Retention", fontsize=18)
axs[0,1].set_xlabel("VDJdb q-score", fontsize=14)
axs[0,1].set_ylabel("Retention", fontsize=14)

# Purity
axs[1,0].bar(x - width/2, height=summ["purity_T"], width=width, tick_label=labels, label="True")
axs[1,0].bar(x + width/2, height=summ["purity_B"], width=width, tick_label=labels, label="Base")
axs[1,0].set_title("Purity", fontsize=18)
axs[1,0].set_xlabel("VDJdb q-score", fontsize=14)
axs[1,0].set_ylabel("Purity", fontsize=14)
axs[1,0].legend()

# Consistency
axs[1,1].bar(x - width/2, height=summ["consistency_T"], width=width, tick_label=labels, label="True")
axs[1,1].bar(x + width/2, height=summ["consistency_B"], width=width, tick_label=labels, label="Base")
axs[1,1].set_title("Consistency", fontsize=18)
axs[1,1].set_xlabel("VDJdb q-score", fontsize=14)
axs[1,1].set_ylabel("Consistency", fontsize=14)
axs[1,1].legend()

fig.suptitle("CDR3 clustering - q-score - VDJdb", fontsize=30)
fig.tight_layout()