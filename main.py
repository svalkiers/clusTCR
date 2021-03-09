from clustcr import Clustering, datasets, read_cdr3, ClusterAnalysis




# cdr3 = read_cdr3('/home/max/Downloads/emerson-2017-natgen/HIP00594.tsv', 'immuneaccess')
cdr3 = datasets.vdjdb_cdr3()

print('1')
Clustering().fit(cdr3)
print('all')
Clustering(n_cpus='all').fit(cdr3)


