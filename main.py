from clustcr import Clustering, datasets, read_cdr3, ClusterAnalysis

cdr3 = datasets.test_cdr3()
features = Clustering().fit(cdr3).compute_features(compute_pgen=True)

ClusterAnalysis(features)


# cdr3 = read_cdr3('/home/max/Downloads/emerson-2017-natgen/HIP00594.tsv', 'immuneaccess')

# print('1')
# Clustering().fit(cdr3)
# print('all')
# Clustering(n_cpus='all').fit(cdr3)


