import pandas as pd
import seaborn as sns
import time
from os.path import join
import matplotlib.pyplot as plt
from clustcr import Clustering, datasets

plt.style.use(['seaborn-white', 'seaborn-paper'])
plt.rc('font', family='serif')
sns.set_palette('Set1')
sns.set_context('paper', font_scale=1.3)


def evaluate_distance_metrics(start, end, step_size, replicates, filename=None):
    final = pd.DataFrame()
    for n in range(start, end, step_size):
        print('###################')
        print(n)
        print('###################')
        for i in range(replicates):
            
            beta = datasets.vdjdb_beta().sample(n)
            epi = datasets.vdjdb_beta(epitopes=True)
            epi = epi[epi.CDR3.isin(beta)]
            
            t = time.time()
            out_HD = Clustering(method='two-step', distance_metric='HAMMING').fit(beta)
            t_hd = time.time() - t
            
            t = time.time()
            out_LD = Clustering(method='two-step', distance_metric='LEVENSHTEIN').fit(beta)
            t_ld = time.time() - t
            
            summ_HD = out_HD.metrics(epi).summary()
            summ_HD['n'] = n
            summ_HD['dm'] = 'Hamming'
            summ_HD['t'] = t_hd
            summ_LD = out_LD.metrics(epi).summary()
            summ_LD['n'] = n
            summ_LD['dm'] = 'Levenshtein'
            summ_LD['t'] = t_ld
            final = final.append(summ_HD)
            final = final.append(summ_LD)
            
    if filename is not None:
        final.to_csv(join('./results/',filename), sep='\t', index=False)
            
    return final
    
def show_results(data):
    colors = sns.color_palette('Set1')
    
    retent = data[data['metrics']=='retention']
    purity = data[data['metrics']=='purity']
    pur_90 = data[data['metrics']=='purity_90']
    consist = data[data['metrics']=='consistency']
    
    fig = plt.figure() 

    gs = fig.add_gridspec(3,2)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[1,0])
    ax4 = fig.add_subplot(gs[1,1])
    ax5 = fig.add_subplot(gs[2,:])
    
    x = data.n.unique()
    ax1.plot(x, retent[retent['dm']=='Hamming'].actual, color=colors[0], label='Hamming')
    ax1.plot(x, retent[retent['dm']=='Levenshtein'].actual, color=colors[1], label='Levenshtein')
    ax1.set_ylabel('Retention')
    
    ax2.plot(x, purity[purity['dm']=='Hamming'].actual, color=colors[0], label='Hamming')
    ax2.plot(x, purity[purity['dm']=='Levenshtein'].actual, color=colors[1], label='Levenshtein')
    ax2.plot(x, purity[purity['dm']=='Hamming'].baseline, color=colors[0], ls='--', lw=.5, label='Hamming - permuted')
    ax2.plot(x, purity[purity['dm']=='Levenshtein'].baseline, color=colors[1], ls='--', lw=.5, label='Levenshtein - permuted')
    ax2.set_ylabel('Purity')
    
    ax3.plot(x, pur_90[pur_90['dm']=='Hamming'].actual, color=colors[0], label='Hamming')
    ax3.plot(x, pur_90[pur_90['dm']=='Levenshtein'].actual, color=colors[1], label='Levenshtein')
    ax3.plot(x, pur_90[pur_90['dm']=='Hamming'].baseline, color=colors[0], ls='--', lw=.5, label='Hamming - permuted')
    ax3.plot(x, pur_90[pur_90['dm']=='Levenshtein'].baseline, color=colors[1], ls='--', lw=.5, label='Levenshtein - permuted')
    ax3.set_ylabel(r'$f_{purity > 0.90}$')
    
    ax4.plot(x, consist[consist['dm']=='Hamming'].actual, color=colors[0], label='Hamming')
    ax4.plot(x, consist[consist['dm']=='Levenshtein'].actual, color=colors[1], label='Levenshtein')
    ax4.plot(x, consist[consist['dm']=='Hamming'].baseline, color=colors[0], ls='--', lw=.5, label='Hamming - permuted')
    ax4.plot(x, consist[consist['dm']=='Levenshtein'].baseline, color=colors[1], ls='--', lw=.5, label='Levenshtein - permuted')
    ax4.set_ylabel('Consistency')
    
    ax5.plot(x, consist[consist['dm']=='Hamming'].t, color=colors[0], label='Hamming')
    ax5.plot(x, consist[consist['dm']=='Levenshtein'].t, color=colors[1], label='Levenshtein')
    ax5.set_ylabel('t (seconds)')
    ax5.set_xlabel('n sequences')
    
    fig.subplots_adjust(top=1.1, hspace=1, wspace=.5)
    
    handles, labels = ax2.get_legend_handles_labels()
    fig.legend(handles, labels, loc='center right',
              bbox_to_anchor=(1.35,0.6))
    
    ax1.text(-0.25, 1.50, 'A', transform=ax1.transAxes,fontsize=20, fontweight='bold', va='top', ha='right')
    ax2.text(-0.25, 1.50, 'B', transform=ax2.transAxes,fontsize=20, fontweight='bold', va='top', ha='right')
    ax3.text(-0.25, 1.50, 'C', transform=ax3.transAxes,fontsize=20, fontweight='bold', va='top', ha='right')
    ax4.text(-0.25, 1.50, 'D', transform=ax4.transAxes,fontsize=20, fontweight='bold', va='top', ha='right')
    ax5.text(-0.1, 1.50, 'E', transform=ax5.transAxes,fontsize=20, fontweight='bold', va='top', ha='right')
     

if __name__=="__main__":
    results = evaluate_distance_metrics(1000,30000,1000,1)
    show_results(results)