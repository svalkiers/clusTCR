import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from clustcr import Clustering, datasets

plt.style.use(['seaborn-white', 'seaborn-paper'])
plt.rc('font', family='serif')
sns.set_palette('Set1')
sns.set_context('paper', font_scale=1.3)


def evaluate_distance_metrics(start, end, step_size, replicates):
    final = pd.DataFrame()
    for n in range(start, end, step_size):
        print('###################')
        print(n)
        print('###################')
        for i in range(replicates):
            beta = datasets.vdjdb_beta().sample(n)
            epi = datasets.vdjdb_beta(epitopes=True)
            epi = epi[epi.CDR3.isin(beta)]
            out_HD = Clustering(method='two-step', distance_metric='HAMMING').fit(beta)
            out_LD = Clustering(method='two-step', distance_metric='LEVENSHTEIN').fit(beta)
            summ_HD = out_HD.metrics(epi).summary()
            summ_HD['n'] = n
            summ_HD['dm'] = 'Hamming'
            summ_LD = out_LD.metrics(epi).summary()
            summ_LD['n'] = n
            summ_LD['dm'] = 'Levenshtein'
            final = final.append(summ_HD)
            final = final.append(summ_LD)
    return final
    
def show_results(data):
    colors = sns.color_palette('Set1')
    metrics = {'retention':'Retention',
               'purity':'Purity',
               'purity_90':r'$f_{purity > 0.90}$',
               'consistency':'Consistency'}
    axes = [(0,0), (0,1), (1,0), (1,1)]
    fig, ax = plt.subplots(2,2)
    hamming = data[data.dm=='Hamming']
    levensh = data[data.dm=='Levenshtein']
    for loc, metric in zip(axes, metrics.keys()):
        hd = hamming[hamming.metrics==metric]
        ld = levensh[levensh.metrics==metric]
        x = hd.n.unique()
        # Hamming y
        y_hd_a = hd.groupby(by='n').mean().actual
        std_hd = hd.groupby(by='n').std().actual
        y_hd_b = hd.groupby(by='n').mean().baseline
        # Levenshtein y
        y_ld_a = ld.groupby(by='n').mean().actual
        std_ld = ld.groupby(by='n').std().actual
        y_ld_b = ld.groupby(by='n').mean().baseline
        ax[loc].plot(x, y_hd_a, label='Hamming', color=colors[0])
        ax[loc].fill_between(x, y_hd_a+std_hd, y_hd_a-std_hd, color=colors[0], label='$\pm$ 1 $\sigma$', alpha=.3)
        ax[loc].plot(x, y_hd_b, label='Baseline', ls='--', color=colors[0])
        ax[loc].plot(x, y_ld_a, label='Levenshtein', color=colors[1])
        ax[loc].fill_between(x, y_ld_a+std_ld, y_ld_a-std_ld, color=colors[1], label='$\pm$ 1 $\sigma$', alpha=.3)
        ax[loc].plot(x, y_ld_b, label='Baseline', ls='--', color=colors[1])
        ax[1,0].set_xlabel("n sequences")
        ax[1,1].set_xlabel("n sequences")
        ax[loc].set_ylabel(metrics[metric])
        
    # plt.subplots_adjust(right=0.75)
    # handles, labels = ax[0,0].get_legend_handles_labels()
    # lgd = fig.legend(handles, labels, loc='center right', bbox_to_anchor=(1.25,0.5))
    
    legend_labels = ['Hamming', 'Baseline', 'Levenshtein',
                     'Baseline', '$\pm$ 1 $\sigma$', '$\pm$ 1 $\sigma$']
    legend = fig.legend([ax[0,0], ax[0,1], ax[1,0], ax[1,1]],
                        labels=legend_labels,
                        loc="lower center",
                        borderaxespad=.1,
                        title="Legend",
                        fontsize="12",
                        ncol=3,
                        bbox_to_anchor=(.5, -0.175)
                       )
    
    fig.tight_layout()
    fig.subplots_adjust(top=1.1)
    ax[0,0].text(-0.225, 1.1, 'A', transform=ax[0,0].transAxes,fontsize=20, fontweight='bold', va='top', ha='right')
    ax[0,1].text(-0.225, 1.1, 'B', transform=ax[0,1].transAxes,fontsize=20, fontweight='bold', va='top', ha='right')
    ax[1,0].text(-0.225, 1.1, 'C', transform=ax[1,0].transAxes,fontsize=20, fontweight='bold', va='top', ha='right')
    ax[1,1].text(-0.225, 1.1, 'D', transform=ax[1,1].transAxes,fontsize=20, fontweight='bold', va='top', ha='right')
    plt.setp(legend.get_title(), fontweight='bold')
    fig.savefig('LD_HD.pdf', format='pdf', bbox_inches='tight', pad_inches=0.3)   

if __name__=="__main__":
    # results = evaluate_distance_metrics(1000,30000,1000,3)
    show_results(results)