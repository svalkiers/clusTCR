from clustcr import datasets, Clustering
import pandas as pd


names = ['HIP00707.tsv', 'HIP00710.tsv', 'HIP13923.tsv', 'HIP14106.tsv']
loaded_data = [datasets.read_cdr3('/home/max/Downloads/emerson-2017-natgen/' + name, 'immuneaccess') for name in names]
concatenated = pd.concat(loaded_data)
print(len(concatenated), 'sequences loaded')


clustering = Clustering(faiss_training_data=loaded_data[0],
                        fitting_data_size=len(concatenated),
                        max_sequence_size=concatenated.str.len().max(),
                        n_cpus='all')
print('trained')

for name, data in zip(names, loaded_data):
    print('precluster', name)
    clustering.batch_precluster(data, name)

print('cluster')
for cluster in clustering.batch_cluster(calc_feature_matrix=True):
    print(cluster.clusters_df)

print(clustering.batch_feature_matrix())

clustering.batch_cleanup()
