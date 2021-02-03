from clustcr import Clustering, datasets
from clustcr.modules.faiss_clustering import datasets as d

POSSIBLE_DATA = [lambda: d.vdj()[0], d.covid19_repertoire]


def max_sequence_size_and_total_fitting_size():
    max_sequence_size = 0
    total_fitting_size = 0
    for data in POSSIBLE_DATA:
        data = data()
        total_fitting_size += len(data)
        data_max_sequence_size = data.str.len().max()
        if data_max_sequence_size > max_sequence_size:
            max_sequence_size = data_max_sequence_size
    return max_sequence_size, total_fitting_size


def batch_cluster():
    max_sequence_size, total_fitting_size = max_sequence_size_and_total_fitting_size()
    cdr3_train = datasets.vdjdb_cdr3().sample(14000)
    clustering = Clustering(faiss_training_data=cdr3_train,
                            max_sequence_size=max_sequence_size,
                            fitting_data_size=total_fitting_size)

    for data in POSSIBLE_DATA:
        cdr3 = data()
        clustering.batch_precluster(cdr3)

    for clusters in clustering.batch_cluster():
        print(clusters.clusters_df)

    clustering.batch_cleanup()


if __name__ == '__main__':
    batch_cluster()
