import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('../results/method_comparison_speed.tsv', sep = '\t')
methods = data.columns[:-1]
n_seq = data[data.columns[-1]]

fig, ax = plt.subplots(figsize=(12,8))
for method in methods:
    ax.scatter(x = n_seq, y = data[method] / 3600, marker = 'D', s = 150, label = method)
    ax.plot(n_seq, data[method] / 3600, alpha = 0.3, lw = 4, ls = '--')

ax.set_title('Speed of CDR3 clustering methods', fontsize = 24)
ax.set_xlabel('Number of input sequences', fontsize = 18)
ax.set_ylabel('Runtime (hours)', fontsize = 18)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)
ax.legend(fontsize = 'xx-large')

fig.savefig('../results/figures/method_comparison_speed.png')