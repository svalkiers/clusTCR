from clustcr import Clustering, metarepertoire
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

datadir = '/home/sebastiaan/PhD/Data/Emerson/emerson_cohort_1/'

plt.style.use(['seaborn-white', 'seaborn-paper'])
plt.rc('font', family='serif')
sns.set_palette('Set1')
sns.set_context('paper', font_scale=1.3)

def evaluate_performance():
    for i in range(1,11):
        sample_size = 10**5 * i
        sample = metarepertoire(datadir,
                                data_format='immuneaccess',
                                n_sequences=sample_size)
        clustering = Clustering(n_cpus=1, method='two-step')
        out = clustering.fit(sample)
    
def show_results():
    colnames = ['method', 'time']
    data = pd.read_csv('/home/sebastiaan/Desktop/time_2.csv', names=colnames)
    sizes = [i for i in range(1,11)] * 3
    data["sizes"] = sizes
    mapping = {"'_faiss'":'Faiss',
               "'MCL_from_preclusters'":'MCL',
               "'_twostep'":'Two-step'}
    data.method = data.method.map(mapping)
    fig, ax = plt.subplots()
    colors = sns.color_palette("Set1")[:2]
    
    x = data.sizes.unique()
    faiss_time = data[data.method=='Faiss'].time.reset_index(drop=True)
    mcl_time = data[data.method=='MCL'].time.reset_index(drop=True)
    total = mcl_time + faiss_time
    
    ax.plot(x, faiss_time, color=colors[0])
    ax.plot(x, total, color=colors[1])
    
    ax.fill_between(x, faiss_time, 0, color=colors[0], alpha=.4, label="Step 1 (Faiss)")
    ax.fill_between(x, total, faiss_time, color=colors[1], alpha=.4, label="Step 2 (MCL)")
    
    # for color, step in zip(colors, data.method.unique()):
    #     x = data.sizes.unique()
    #     y = data[data.method==step].time
    #     ax.plot(x, y, label=step, color=color)
    #     ax.scatter(x, y)
    ax.legend()
    ax.set_xlabel(r"n sequences ($\times 10^5$)")
    ax.set_ylabel("Runtime (s)")
    ax.set_xlim(1,10)
    ax.set_ylim(0,1000)
    fig.tight_layout()
    fig.savefig("./results/figures/step_comparison.svg", format='svg')
    
    return faiss_time, mcl_time
    

a, b = show_results()