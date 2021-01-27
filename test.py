from src import Clustering, datasets

test_cdr3 = datasets.test_cdr3()
clustering = Clustering(method='faiss').fit(test_cdr3)
print(clustering.clusters_df)


