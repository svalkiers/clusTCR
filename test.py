from clustcr import Clustering, datasets

cdr3 = datasets.test_cdr3()
clustering = Clustering().fit(cdr3)
print(clustering.clusters_df)




