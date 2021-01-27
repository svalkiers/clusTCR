from src import Clustering, datasets

print(Clustering(method='faiss').fit(datasets.test_cdr3()).clusters_df)
