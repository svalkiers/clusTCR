import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import networkx as nx
import os

from clustcr.clustering.clustering import ClusteringResult
from clustcr.input.vdjdb import parse_vdjdb
from clustcr.modules.gliph2.gliph2 import GLIPH2
from clustcr.modules.ismart.ismart import iSMART

from tcrdist.repertoire import TCRrep
from tcrdist.rep_funcs import  compute_pw_sparse_out_of_memory
from sklearn.cluster import DBSCAN, KMeans
from scipy.spatial.distance import cdist


def TCRDist(data=None, chain='beta', sparse=False):
    
    if data is None:    
        vdjdb = parse_vdjdb(os.path.abspath('./clustcr/input/vdjdb/vdjdb_full.txt'), q=1)
    else:
        vdjdb = data

    if chain == 'beta':
        cdr3 = 'cdr3_b_aa'
        v_name = 'v_b_gene'
        vdjdb = vdjdb.drop(columns=['cdr3.alpha', 'v.alpha'])
        vdjdb = vdjdb.rename(columns={'cdr3.beta':cdr3,
                                      'v.beta':v_name})
    elif chain == 'alpha':
        cdr3 = 'cdr3_a_aa'
        v_name = 'v_a_gene'
        vdjdb = vdjdb.drop(columns=['cdr3.beta', 'v.beta'])
        vdjdb = vdjdb.rename(columns={'cdr3.alpha':cdr3,
                                      'v.alpha':v_name})        
    
    df_epi = vdjdb[[cdr3, v_name, 'antigen.epitope']].dropna().drop_duplicates()
    seq = df_epi.drop(columns = ['antigen.epitope']).drop_duplicates().reset_index(drop=True)
    
    # gt = df_epi.drop(columns = [v_name])
    gt = df_epi.rename(columns = {cdr3:'CDR3',
                                  v_name:'V',
                                  'antigen.epitope':'Epitope'})
    
    final = pd.DataFrame()
    
    if sparse:
        
        # tr = TCRrep(cell_df = seq,
        #     organism = 'human',
        #     chains = ['beta'],
        #     db_file = 'alphabeta_gammadelta_db.tsv',
        #     compute_distances = False,
        #     store_all_cdr = False)

        # S, fragments = compute_pw_sparse_out_of_memory( tr = tr,
        #     row_size      = 500,
        #     pm_processes  = 8,
        #     pm_pbar       = True,
        #     max_distance  = 500,
        #     matrix_name   = 'rw_beta',
        #     reassemble    = True,
        #     cleanup       = True)
        
        tr = TCRrep(cell_df = seq,
                    organism = 'human',
                    chains = ['beta'],
                    db_file = 'alphabeta_gammadelta_db.tsv',
                    compute_distances = False)
        
        tr.cpus = 2
        tr.compute_sparse_rect_distances(radius = 200, chunk_size = 500)
        S = tr.rw_beta
        
    else:
        
        tr = TCRrep(cell_df = seq,
            organism = 'human',
            chains = [chain],
            db_file = 'alphabeta_gammadelta_db.tsv',
            compute_distances = True)
    
        S = tr.pw_cdr3_b_aa
        
    return S, seq, gt

def normalize(edge):
    n1, n2 = edge
    if n1 > n2:
        n1, n2 = n2, n1
    return (n1, n2)

def greedy_clustering(dm, threshold):

    edges = np.argwhere(dm<=threshold)
    print(len(edges))
    edges = set(map(normalize, edges)) # Remove duplicated edges
    edges = np.array(list(edges)) # Edgelist to array
    print(len(edges))
    
    cid = 0
    res = pd.DataFrame()
    
    while len(edges) > 0:
        
        # print(len(edges))
        
        G = nx.from_edgelist(edges)
        degrees = pd.DataFrame(G.degree(), columns=['node', 'degree'])
        degrees = degrees.set_index('node')
        degrees = degrees.sort_values(by='degree', ascending=False)
        max_degree = degrees.idxmax().values
        
        cluster = edges[np.where(edges[:,0]==max_degree)]
        ids = np.unique(cluster)
        cids = [cid] * len(ids)

        if len(ids) <= 1:
            break
        
        cluster_iter = pd.DataFrame({'seq_id':ids,'cluster':cids})
        res = res.append(cluster_iter)
        
        for i in ids:
            edges = edges[np.where(edges[:,0]!=i)] # Remove from column 1
            edges = edges[np.where(edges[:,1]!=i)] # Remove from column 2
        
        cid += 1
            
    return res

def cluster_TCRDist_matrix(dm=None, cdr3=None, gt=None, method='DBSCAN'):
    
    methods = ['GREEDY', 'DBSCAN', 'KMEANS']
    assert method.upper() in methods, r'Please choose one of the following: /n %s' % methods
    
    if dm is None:
        dm, cdr3, gt = TCRDist(sparse=False)
    if gt is None:
        dm, cdr3, gt = TCRDist(sparse=False)
    if cdr3 is None:
        dm, cdr3, gt = TCRDist(sparse=False)
        
    
        
    if method.upper() == 'GREEDY':
        
        # Greedy clustering
        clusters = greedy_clustering(dm, 12)
        clusters = clusters.rename(columns={'seq_id':'Index'})
        clusters = clusters.set_index('Index', drop=True)
        clusters = clusters.merge(right=cdr3, left_index=True, right_index=True)
        clusters = clusters.rename(columns={'cdr3_b_aa':'CDR3',
                                            'v_b_gene':'V'})

        metrics = ClusteringResult(clusters).metrics(gt)
        return metrics.summary()
    
    elif method.upper() == 'DBSCAN':
        
        # DBSCAN
        clustering = DBSCAN(eps=250, min_samples=2, n_jobs=-1).fit(dm)
        labels = clustering.labels_
        clusters = cdr3.rename(columns={'cdr3_b_aa':'CDR3',
                                        'v_b_gene':'V'})
        clusters['cluster'] = labels
        clusters = clusters[clusters['cluster']!=-1]
        metrics = ClusteringResult(clusters).metrics(gt)
        return metrics.summary()
    
    else:
        
        # K-Means
        kmeans = KMeans(n_clusters=500).fit(dm)
        labels = kmeans.labels_
        clusters = cdr3.rename(columns={'cdr3_b_aa':'CDR3',
                                        'v_b_gene':'V'})
        clusters['cluster'] = labels
        metrics = ClusteringResult(clusters).metrics(gt)
        return metrics.summary()

# d_reg, s, g = TCRDist()
d_spa, s, g = TCRDist(sparse=True)       
# res = cluster_TCRDist_matrix(d_spa.todense(), s, g)


# # KMeans clustering
# final_summ = pd.DataFrame()
# for k in range(50,501,50):
#     print(k)
#     kmeans = KMeans(n_clusters=k).fit(dm)
#     labels = kmeans.labels_
#     clusters = cdr3.rename(columns={'cdr3_b_aa':'CDR3',
#                                     'v_b_gene':'V'})
#     clusters['cluster'] = labels
#     metrics = ClusteringResult(clusters).metrics(gt)
#     summ = metrics.summary()
#     summ['k'] = k
#     final_summ = final_summ.append(summ)

# ssd = []
# n_clusters = []
# K = range(5,301,5)
# for k in K:
#     print(k)
#     kmeans = KMeans(n_clusters=k).fit(dm)
#     ssd.append(kmeans.inertia_)
#     n_clusters.append(k) 
    
#     ax.plot(n_clusters, ssd)
#     fig.tight_layout()
#     plt.show()
    
# fig, ax = plt.subplots()
# ax.plot(n_clusters, ssd, 'bx-')
# ax.set_xlabel('k')
# ax.set_ylabel('Sum of squared distances')
# fig.tight_layout()
# fig.savefig('./results/figures/tcrdist_kmeans_elbow.eps', format='eps')


# plt.style.use(['seaborn-white', 'seaborn-paper'])
# plt.rc('font', family='serif')
# sns.set_palette('Set1')
# sns.set_context('paper', font_scale=1.3)

# clrs = sns.color_palette("Set1")
# fig, ax = plt.subplots(nrows=2, ncols=2)

# thresholds = final.k.unique()
# # retention = final[final.metrics=='retention'].actual
# purity_a = final[final.metrics=='purity'].actual
# purity_b = final[final.metrics=='purity'].baseline
# pur_90_a = final[final.metrics=='purity_90'].actual
# pur_90_b = final[final.metrics=='purity_90'].baseline
# consistency_a = final[final.metrics=='consistency'].actual
# consistency_b = final[final.metrics=='consistency'].baseline

# # plt.setp(ax, xlim=(1,30), ylim=(0,1), xlabel='k')

# ax[0,0].plot(n_clusters, ssd, color=clrs[0])
# ax[0,0].scatter(n_clusters, ssd, color=clrs[0], marker='x')
# ax[0,0].set_ylabel('Sum of squared distances')

# ax[0,1].plot(thresholds, purity_a, color=clrs[0], label='Actual')
# ax[0,1].plot(thresholds, purity_b, color=clrs[1], label='Baseline', ls='--')
# ax[0,1].set_ylabel('Purity')
# ax[0,1].legend()

# ax[1,0].plot(thresholds, pur_90_a, color=clrs[0], label='Actual')
# ax[1,0].plot(thresholds, pur_90_b, color=clrs[1], label='Baseline', ls='--')
# ax[1,0].set_ylabel(r'$f_{purity > 0.90}$')
# ax[1,0].legend()

# ax[1,1].plot(thresholds, consistency_a, color=clrs[0], label='Actual')
# ax[1,1].plot(thresholds, consistency_b, color=clrs[1], label='Baseline', ls='--')
# ax[1,1].set_ylabel('Consistency')
# ax[1,1].legend()

# ax[0,0].text(-0.175, 1.1, 'A', transform=ax[0,0].transAxes,fontsize=20, fontweight='bold', va='top', ha='right')
# ax[0,1].text(-0.175, 1.1, 'B', transform=ax[0,1].transAxes,fontsize=20, fontweight='bold', va='top', ha='right')
# ax[1,0].text(-0.175, 1.1, 'C', transform=ax[1,0].transAxes,fontsize=20, fontweight='bold', va='top', ha='right')
# ax[1,1].text(-0.175, 1.1, 'D', transform=ax[1,1].transAxes,fontsize=20, fontweight='bold', va='top', ha='right')

# # handles, labels = ax[0,1].get_legend_handles_labels()
# # fig.legend(handles, labels, loc='center right', bbox_to_anchor=(1.25, .8))

# fig.subplots_adjust(right=1, hspace=0.25, wspace=0.35, top=1.5)
# fig.savefig('./results/figures/tcrdist_kmeans.eps', format='eps', bbox_inches='tight')
# plt.show()

# clusters = cdr3.rename(columns={'cdr3_b_aa':'CDR3',
#                                 'v_b_gene':'V'})
# clusters['cluster'] = labels
# clusters = clusters[clusters['cluster']!=-1]
# metrics = ClusteringResult(clusters).metrics(gt)
# print(metrics.summary())



def plot_results(final, figname='tcrdist_dbscan'):
          
    plt.style.use(['seaborn-white', 'seaborn-paper'])
    plt.rc('font', family='serif')
    sns.set_palette('Set1')
    sns.set_context('paper', font_scale=1.3)
            
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, sharex=True)
    
    colors = sns.color_palette('Set1')
    
    purity = final[final.metrics=='purity']
    purity_90 = final[final.metrics=='purity_90']
    consistency = final[final.metrics=='consistency']
    retention = final[final.metrics=='retention'].actual
    
    ax1.plot(retention, purity.actual, c=colors[0], label='Actual')
    ax1.plot(retention, purity.baseline, c= colors[1], label='Baseline', ls='--')
    ax1.set_ylabel('Purity')
    ax1.set_xlim(retention.min(), retention.max())
    
    ax2.plot(retention, purity_90.actual, c=colors[0], label='Actual')
    ax2.plot(retention, purity_90.baseline, c= colors[1], label='Baseline', ls='--')
    ax2.set_ylabel(r'$f_{purity > 0.90}$')
    
    ax3.plot(retention, consistency.actual, c=colors[0], label='Actual')
    ax3.plot(retention, consistency.baseline, c= colors[1], label='Baseline', ls='--')
    ax3.set_xlabel('Retention')
    ax3.set_ylabel('Consistency')
    
    ax1.text(-0.1, 1.1, 'A', transform=ax1.transAxes,fontsize=20, fontweight='bold', va='top', ha='right')
    ax2.text(-0.1, 1.1, 'B', transform=ax2.transAxes,fontsize=20, fontweight='bold', va='top', ha='right')
    ax3.text(-0.1, 1.1, 'C', transform=ax3.transAxes,fontsize=20, fontweight='bold', va='top', ha='right')
    
    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, loc='center right', bbox_to_anchor=(1.25, .8))
    
    fig.subplots_adjust(right=1, hspace=0.25, top=1.5)
    fig.savefig(figname+'.eps', format='eps', bbox_inches='tight')
    # fig.tight_layout()
    


# plot_results(final)


# fig, ax = plt.subplots()

# colors = sns.color_palette('Set1')

# purity = final[final.metrics=='purity']
# purity_90 = final[final.metrics=='purity_90']
# consistency = final[final.metrics=='consistency']
# retention = final[final.metrics=='retention'].actual

# ax.plot(retention, purity.actual, c=colors[0], label='Actual')
# ax.plot(retention, purity.baseline, c= colors[0], label='Baseline', ls='--')
# # ax.set_ylabel('Purity')
# ax.set_xlim(retention.min(), retention.max())

# ax.plot(retention, purity_90.actual, c=colors[1], label='Actual')
# ax.plot(retention, purity_90.baseline, c= colors[1], label='Baseline', ls='--')
# # ax.set_ylabel(r'$f_{purity > 0.90}$')

# ax.plot(retention, consistency.actual, c=colors[2], label='Actual')
# ax.plot(retention, consistency.baseline, c= colors[2], label='Baseline', ls='--')
# ax.set_xlabel('Retention')
# # ax.set_ylabel('Consistency')

# # handles, labels = ax1.get_legend_handles_labels()
# # fig.legend(handles, labels, loc='center right', bbox_to_anchor=(1.25, .5))

# fig.subplots_adjust(right=1.)
# fig.tight_layout()
