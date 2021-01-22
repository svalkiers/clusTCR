import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('../results/method_comparison_accuracy.tsv', sep = '\t')

metrics = ['Retention', 'Purity', 'Consistency', 'avg_cluster_size']
for metric in metrics:
    fig, ax = plt.subplots(figsize=(12,8))
    for method in df.Method.unique():
        method_df = df[df['Method'] == method]
        x = method_df.Q
        y = method_df[metric]
        ax.scatter(x, y, marker = 'D', s = 150, label = method)
        ax.plot(x, y, alpha = 0.3, lw = 4, ls = '--')
    
    ax.set_title('{} of TCR clustering approaches'.format(metric), fontsize = 24)
    ax.set_xlabel('VDJdb subset', fontsize = 18, labelpad=20)
    ax.set_ylabel('{}'.format(metric), fontsize = 18)
    ax.tick_params(axis='y', which='major', labelsize=12)
    ax.tick_params(axis='x', which='major', labelsize=18)
    # ax.tick_params(axis='both', which='minor', labelsize=8)
    ax.legend(fontsize = 'xx-large')
    
    locs = [0, 1, 2]
    labels = ['All', 'Q $\geq$ 1', 'Q $\geq$ 2']
    plt.xticks(ticks=locs, labels=labels)
    plt.grid(axis='y')