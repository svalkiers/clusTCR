import pandas as pd
import sys
import os

path = os.getcwd() + '/load_files/'
sys.path.insert(1, path)

from datasets import vdj_cdr3_small, vdj_gliph2_small, vdj_epitopes_small
from gliph2 import GLIPH2
from ismart import iSMART
from clusTCR import Clustering, Metrics
from load_files.import_functions import path_in_data

summary = {'Method' : [], 'Q' : [], 'Retention' : [], 'Purity' : [], 'Consistency' : []}

# Different quality subsets of VDJdb
qscore = [0, 1, 2]
for q in qscore:
    
    # Data
    epitopes = vdj_epitopes_small()
    data_gliph2 = vdj_gliph2_small(q = q)
    data_cdr3 = vdj_cdr3_small(q = q)
    
    # GLIPH2
    summary['Method'].append('GLIPH2')
    summary['Q'].append(q)
    GLIPH2(data_gliph2)
    gliph2_output = pd.read_csv(path_in_data('methods/GLIPH2_{}.tsv'.format(len(data_gliph2))))
    gliph2_metrics = Metrics(gliph2_output, epitopes)
    summary['Retention'].append(gliph2_metrics.retention())
    summary['Purity'].append(gliph2_metrics.purity())
    summary['Consistency'].append(gliph2_metrics.consistency())
    
    # iSMART
    summary['Method'].append('iSMART')
    summary['Q'].append(q)
    iSMART(data_cdr3)
    ismart_output = pd.read_csv(path_in_data('methods/iSMART_{}.tsv'.format(len(data_cdr3))))
    ismart_metrics = Metrics(ismart_output, epitopes)
    summary['Retention'].append(ismart_metrics.retention())
    summary['Purity'].append(ismart_metrics.purity())
    summary['Consistency'].append(ismart_metrics.consistency())
    
    # Two-step
    summary['Method'].append('Two-step')
    summary['Q'].append(q)
    twostep = Clustering(data_cdr3, method = 'two-step')
    ts_output = twostep.cluster_sequences()
    ts_metrics = Metrics(ts_output, epitopes)
    summary['Retention'].append(ts_metrics.retention())
    summary['Purity'].append(ts_metrics.purity())
    summary['Consistency'].append(ts_metrics.consistency())

summary = pd.DataFrame(summary)    
summary.to_csv('../results/method_comparison_accuracy.tsv', sep = '\t', index = False)