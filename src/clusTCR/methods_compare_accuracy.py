import pandas as pd

from load_files.datasets import vdj_cdr3_small, vdj_gliph2_small, vdj_epitopes_small, vdj_tcrdist_small
from gliph2 import GLIPH2
from ismart import iSMART
from clusTCR import Clustering, Metrics
from load_files.import_functions import path_in_data
from tcrdist.repertoire import TCRrep
from sklearn.cluster import DBSCAN

summary = {'Method' : [], 
           'Q' : [], 
           'Retention' : [], 
           'Purity' : [], 
           'Consistency' : [], 
           'avg_cluster_size' : []}

# Different quality subsets of VDJdb
qscore = [0, 1, 2]
for q in qscore:
    
    # Data
    epitopes = vdj_epitopes_small(q = q)
    data_gliph2 = vdj_gliph2_small(q = q)
    data_cdr3 = vdj_cdr3_small(q = q)
    data_td = vdj_tcrdist_small(q = q)
    
    # GLIPH2
    summary['Method'].append('GLIPH2')
    summary['Q'].append(q)
    GLIPH2(data_gliph2)
    gliph2_output = pd.read_csv(path_in_data('methods/GLIPH2_{}.tsv'.format(len(data_gliph2))), sep='\t')
    gliph2_metrics = Metrics(gliph2_output, epitopes)
    summary['Retention'].append(gliph2_metrics.retention())
    summary['Purity'].append(gliph2_metrics.purity()['Actual'])
    summary['Consistency'].append(gliph2_metrics.consistency()['Actual'])
    summary['avg_cluster_size'].append(gliph2_output.cluster.value_counts().mean())
    
    # iSMART
    summary['Method'].append('iSMART')
    summary['Q'].append(q)
    iSMART(data_cdr3)
    ismart_output = pd.read_csv(path_in_data('methods/iSMART_{}.tsv'.format(len(data_cdr3))), sep='\t')
    ismart_metrics = Metrics(ismart_output, epitopes)
    summary['Retention'].append(ismart_metrics.retention())
    summary['Purity'].append(ismart_metrics.purity()['Actual'])
    summary['Consistency'].append(ismart_metrics.consistency()['Actual'])
    summary['avg_cluster_size'].append(ismart_output.cluster.value_counts().mean())
    
    # Two-step
    summary['Method'].append('Two-step')
    summary['Q'].append(q)
    twostep = Clustering(data_cdr3, method = 'two-step', n_cpus=8)
    ts_output = twostep.cluster_sequences()
    ts_metrics = Metrics(ts_output, epitopes)
    summary['Retention'].append(ts_metrics.retention())
    summary['Purity'].append(ts_metrics.purity()['Actual'])
    summary['Consistency'].append(ts_metrics.consistency()['Actual'])
    summary['avg_cluster_size'].append(ts_output.cluster.value_counts().mean())
    
    # tcrdist3
    if q != 0:
        print('pw dist calculations for {} sequences with tcrdist.'.format(len(data_td)))
        summary['Method'].append('tcrdist3')
        summary['Q'].append(q)
        tr = TCRrep(cell_df = data_td,
                    organism = 'human',
                    chains = ['beta'],
                    db_file='alphabeta_gammadelta_db.tsv')
        d = tr.pw_cdr3_b_aa
        clustering = DBSCAN(eps=18, min_samples=2, n_jobs=-1).fit(d)
        labels = clustering.labels_
        data_td['cluster'] = labels
        data_td = data_td[data_td['cluster']!=-1]
        data_td.rename(columns={'cdr3_b_aa':'CDR3','v_b_gene':'V'}, inplace=True)
        td_metrics = Metrics(data_td, epitopes)
        summary['Retention'].append(td_metrics.retention())
        summary['Purity'].append(td_metrics.purity()['Actual'])
        summary['Consistency'].append(td_metrics.consistency()['Actual'])
        summary['avg_cluster_size'].append(data_td.cluster.value_counts().mean())
    else:
        summary['Method'].append('tcrdist3')
        summary['Q'].append(q)
        summary['Retention'].append(None)
        summary['Purity'].append(None)
        summary['Consistency'].append(None)
        summary['avg_cluster_size'].append(None)

summary = pd.DataFrame(summary)    
summary.to_csv('../results/method_comparison_accuracy.tsv', sep = '\t', index = False)