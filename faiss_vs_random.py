import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from clustcr import Clustering, datasets

plt.style.use(['seaborn-white', 'seaborn-paper'])
plt.rc('font', family='serif')
sns.set_palette('Set1')
sns.set_context('paper', font_scale=1.3)

beta = datasets.vdjdb_beta()
epitopes = datasets.vdjdb_beta(epitopes=True)

def evaluate_faiss(data, start, end, step_size):
    final = pd.DataFrame()
    for s in range(start,end,step_size):
        out = Clustering(method='faiss', faiss_cluster_size=s).fit(data)
        metrics = out.metrics(epitopes)
        summ = metrics.summary()
        summ["s_cluster"] = s
        final = final.append(summ)
    return final

def evaluate_random(data, start, end, step_size):
    final = pd.DataFrame()
    for s in range(start,end,step_size):
        out = Clustering(method='random', rnd_chunk_size=s).fit(data)
        metrics = out.metrics(epitopes)
        summ = metrics.summary()
        summ["s_cluster"] = s
        final = final.append(summ)
    return final
    
# final.to_csv('faiss_evaluation_2_100.csv', index=False)

final_faiss = evaluate_faiss(beta, 100, 10001, 100)
final_random = evaluate_random(beta, 100, 10001, 100)

x = final_faiss.s_cluster.unique()
purity_a = final_faiss[final_faiss.metrics=='purity'].actual
purity_r = final_random[final_random.metrics=='purity'].actual
p90_a = final_faiss[final_faiss.metrics=='purity_90'].actual
p90_r = final_random[final_random.metrics=='purity_90'].actual
cons_a = final_faiss[final_faiss.metrics=='consistency'].actual
cons_r = final_random[final_random.metrics=='consistency'].actual

fig, (ax1, ax3) = plt.subplots(nrows=2,ncols=1,sharex=True)

ax1.plot(x, purity_a, label='Faiss')
ax1.plot(x, purity_r, ls='--', label='Random')
ax1.set_ylabel('Purity')
# ax1.legend()

# ax2.plot(x, p90_a, label='faiss')
# ax2.plot(x, p90_r, ls='--', label='random')
# ax2.set_ylabel(r'$f_{p\geq0.90}$')
# ax2.legend()

ax3.plot(x, cons_a, label='faiss')
ax3.plot(x, cons_r, ls='--', label='random')
ax3.set_xlabel('Average cluster size')
ax3.set_ylabel('Consistency')
# ax3.legend()

ax1.text(-0.1, 1.1, 'A', transform=ax1.transAxes,fontsize=20, fontweight='bold', va='top', ha='right')
ax3.text(-0.1, 1.1, 'B', transform=ax3.transAxes,fontsize=20, fontweight='bold', va='top', ha='right')

handles, labels = ax1.get_legend_handles_labels()
fig.legend(handles, labels, loc='center right')
plt.subplots_adjust(right=0.80)
fig.savefig("./results/figures/faiss_vs_random.eps", format="eps")
# fig.tight_layout()