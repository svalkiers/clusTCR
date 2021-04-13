import pandas as pd
import matplotlib.pyplot as plt

colnames = ['method', 'time']
t = pd.read_csv('/home/sebastiaan/Desktop/time.csv', names=colnames)

s = []
for i in range(1000,10001,1000):
    s += [i] * 3
c1 = [1] * len(s)
c2 = [2] * len(s)
c3 = [3] * len(s)
c = c1 + c2 + c3
s = s * 3
t['s'] = s
t['c'] = c

mapping = {"'hashing_method'":'hashing',
           "'hamming_method'":'hamming-alt',
           "'levenshtein_method'":'levenshtein'}

for d in t.c.unique():
    dist = t[t.c==d]
    dist.method = dist.method.map(mapping)
    fig, ax = plt.subplots()
    x = dist.s.unique()
    for meth in dist.method.unique():
        ax.plot(x, dist[dist.method==meth].time, label=meth)
    ax.set_xlabel('n sequences')
    ax.set_ylabel('time (seconds)')
    ax.set_title(f'Cut-off = {d}')
    ax.legend()
    fig.tight_layout()
