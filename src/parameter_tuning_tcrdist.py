import pandas as pd

from tcrdist.repertoire import TCRrep
from load_files.datasets import vdj_gliph2_small, vdj_epitopes_small
from matplotlib import pyplot as plt
from sklearn.cluster import DBSCAN
from clusTCR import Metrics

# Prepare data
vdj = vdj_gliph2_small()
vdj.drop(columns = ['subject', 'count'], inplace = True)
vdj.drop_duplicates(inplace = True)
tcr = vdj.rename(columns = {'CDR3' : 'cdr3_b_aa', 'V' : 'v_b_gene'})

# Calculate pw distances
tr = TCRrep(cell_df = tcr,
            organism = 'human',
            chains = ['beta'],
            db_file='alphabeta_gammadelta_db.tsv')

d = tr.pw_cdr3_b_aa

# Clustering
res = []
for dist in [*range(10,500,10)]:
    print('DBSCAN clustering with eps: %s' % dist)
    clustering = DBSCAN(eps=dist, min_samples=2, n_jobs=-1).fit(d)
    labels = clustering.labels_
    
    # Evaluation of clustering results
    vdj['cluster'] = labels
    vdj_copy = vdj.copy()
    vdj_copy = vdj_copy[vdj_copy['cluster']!=-1]
    
    epitopes = vdj_epitopes_small()
    metrics = Metrics(vdj_copy, epitopes)
    res.append(metrics.summary())

stats = pd.DataFrame(res, columns = ['retention', 'purity', 'consistency'])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,8))
ax1.plot(stats.retention, stats.purity)
ax2.plot(stats.retention, stats.consistency)

ax1.set_title('', fontsize = 24)
ax1.set_xlabel('Retention', fontsize = 18)
ax1.set_ylabel('Purity', fontsize = 18)
ax1.tick_params(axis='both', which='major', labelsize=12)
ax1.tick_params(axis='both', which='minor', labelsize=8)

ax2.set_title('', fontsize = 24)
ax2.set_xlabel('Retention', fontsize = 18)
ax2.set_ylabel('Consistency', fontsize = 18)
ax2.tick_params(axis='both', which='major', labelsize=12)
ax2.tick_params(axis='both', which='minor', labelsize=8)

plt.suptitle('DBSCAN results tcrdist matrix', fontsize=30)