# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 11:32:53 2020

@author: Sebastiaan Valkiers

Contact: sebastiaan.valkiers@uantwerpen.be
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

metrics = ["Purity", "Consistency"]
datasets = ["VDJdb", "Dash"]

for dataset in datasets:
    
    # Infile
    infile = "../results/bm_metrics_" + dataset + ".tsv"
    
    # Read infile and create df per Expansion param value
    df = pd.read_csv(infile, sep="\t")
    df_2 = df[df["Expansion"]==2]
    df_3 = df[df["Expansion"]==3]
    df_4 = df[df["Expansion"]==4]
    df_5 = df[df["Expansion"]==5]

    for metric in metrics:
    
        fig, ax = plt.subplots(figsize=[12,8])
        
        # Regular
        ax.plot(df_2["Inflation"], df_2[metric + "_r"], c="b", label="Expansion = 2")
        ax.plot(df_3["Inflation"], df_3[metric + "_r"], c="r", label="Expansion = 3")
        ax.plot(df_4["Inflation"], df_4[metric + "_r"], c="g", label="Expansion = 4")
        ax.plot(df_5["Inflation"], df_5[metric + "_r"], c="y", label="Expansion = 5")
        
        # Baseline
        ax.plot(df_2["Inflation"], df_2[metric + "_b"], "--", c="b", label="Baseline, Expansion = 2")
        ax.plot(df_3["Inflation"], df_3[metric + "_b"], "--", c="r", label="Baseline, Expansion = 3")
        ax.plot(df_4["Inflation"], df_4[metric + "_b"], "--", c="g", label="Baseline, Expansion = 4")
        ax.plot(df_5["Inflation"], df_5[metric + "_b"], "--", c="y", label="Baseline, Expansion = 5")
        
        # Aesthetics
        ax.set_ylabel(metric, fontsize=16)
        ax.set_xlabel("Inflation", fontsize=16)
        ax.set_title("MCL hyperparameters versus " + metric.lower() + " - " + dataset, fontsize=18)
        ax.legend()
        

# Summary plot comparing different q-score cutoffs (VDJdb)
summ = pd.read_csv("../results/summary_VDJdb.tsv", sep="\t")

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
axs[1,0].bar(x - width/2, height=summ["purity_r"], width=width, tick_label=labels, label="True")
axs[1,0].bar(x + width/2, height=summ["purity_b"], width=width, tick_label=labels, label="Baseline")
axs[1,0].set_title("Purity", fontsize=18)
axs[1,0].set_xlabel("VDJdb q-score", fontsize=14)
axs[1,0].set_ylabel("Purity", fontsize=14)
axs[1,0].legend()

# Consistency
axs[1,1].bar(x - width/2, height=summ["consistency_r"], width=width, tick_label=labels, label="True")
axs[1,1].bar(x + width/2, height=summ["consistency_b"], width=width, tick_label=labels, label="Baseline")
axs[1,1].set_title("Consistency", fontsize=18)
axs[1,1].set_xlabel("VDJdb q-score", fontsize=14)
axs[1,1].set_ylabel("Consistency", fontsize=14)
axs[1,1].legend()

fig.suptitle("CDR3 clustering summary - VDJdb", fontsize=30)
fig.tight_layout()