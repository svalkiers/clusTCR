import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('../results/methods_benchmarked.tsv', sep = '\t')

fig, ax = plt.subplots()
for method in df.method.unique():
    method_df = df[df['method'] == method]
    x = method_df.n_sequences
    y = method_df.runtime
    ax.scatter(x, y, label = method)
    ax.plot(x, y, alpha = 0.3)

ax.set_xlabel('Number of input sequences')
ax.set_ylabel('Runtime (seconds)')
ax.set_title('TCR clustering methods - speed comparison')
ax.legend()

fig, ax = plt.subplots()
for method in df.method.unique():
    method_df = df[df['method'] == method]
    x = method_df.n_sequences
    y = method_df.purity
    ax.scatter(x, y, label = method)
    ax.plot(x, y, alpha = 0.3)

ax.set_xlabel('Number of input sequences')
ax.set_ylabel('Purity')
ax.set_title('TCR clustering methods - accuracy comparison (1)')
ax.legend()

fig, ax = plt.subplots()
for method in df.method.unique():
    method_df = df[df['method'] == method]
    x = method_df.n_sequences
    y = method_df.consistency
    ax.scatter(x, y, label = method)
    ax.plot(x, y, alpha = 0.3)

ax.set_xlabel('Number of input sequences')
ax.set_ylabel('Consistency')
ax.set_title('TCR clustering methods - accuracy comparison (2)')
ax.legend()

fig, ax = plt.subplots()
for method in df.method.unique():
    method_df = df[df['method'] == method]
    x = method_df.n_sequences
    y = method_df.retention
    ax.scatter(x, y, label = method)
    ax.plot(x, y, alpha = 0.3)

ax.set_xlabel('Number of input sequences')
ax.set_ylabel('Retention')
ax.set_title('TCR clustering methods - accuracy comparison (3)')
ax.legend()