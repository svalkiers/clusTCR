import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data = pd.read_csv('../results/method_comparison_speed.tsv', sep = '\t')
methods = data.columns[:-1]
n_seq = data[data.columns[-1]]

fig, ax = plt.subplots(figsize=(12,8))
for method in methods:
    ax.scatter(n_seq, y = np.log10(data[method]), marker = 'D', s = 150, label = method)
    ax.plot(n_seq, np.log10(data[method]), alpha = 0.3, lw = 4, ls = '--')

ax.set_title('Speed of CDR3 clustering methods', fontsize = 24)
ax.set_xlabel('Number of input sequences', fontsize = 18)
ax.set_ylabel('Runtime (log-scaled)', fontsize = 18)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)
ax.legend(fontsize = 'xx-large')
# locs, labels = plt.yticks()
locs = np.log10([1,60,3600,86400])
labels = ["1 sec", "1 min", "1 hr", "1 day"]
plt.ylim((-0.5,5.5))
plt.yticks(ticks=locs, labels=labels)
plt.grid(axis='y')

fig.savefig('../results/figures/method_comparison_speed_log.png')