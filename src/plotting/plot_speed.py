import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('../results/methods_compared_speed.tsv', sep = '\t')
methods = data.columns[:-1]
n_seq = data[data.columns[-1]]

fig, ax = plt.subplots()
for method in methods:
    ax.scatter(x = n_seq, y = data[method], label = method)
    ax.plot(x = n_seq, y = data[method])

ax.set_title('Speed of CDR3 clustering methods')
ax.set_xlabel('Number of input sequences')
ax.set_ylabel('Runtime (seconds)')
ax.legend()