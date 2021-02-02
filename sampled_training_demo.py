from clustcr import Clustering, datasets
import time

cdr3 = datasets.vdjdb_cdr3()


def sampled_training():
    print('Sampled training')
    t = time.time()
    max_sequence_size = cdr3.str.len().max()
    cdr3_train = cdr3.sample(14000)
    clustering = Clustering(faiss_training_data=cdr3_train,
                            max_sequence_size=max_sequence_size,
                            fitting_data_size=len(cdr3))
    result = clustering.fit(cdr3)
    print(result.clusters_df)
    print('in', time.time() - t)


def full_training():
    print('Full training')
    t = time.time()
    clustering = Clustering()
    result = clustering.fit(cdr3)
    print(result.clusters_df)
    print('in', time.time() - t)


sampled_training()
full_training()
