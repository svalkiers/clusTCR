import multiprocessing
import parmap
import pandas as pd
import time
from collections import defaultdict


class GreedyClustering:
    """
    Performs greedy clustering on a Hamming Distance of 1.
    Meaning, when a network of a set of CDR3 sequences is constructed with an edge for each pair that
    has a Hamming Distance of 1, each connected component in the graph results in a cluster.
    These edges can either be constructed using an
    optimized hashing (use_hashing=True) or parallelized (use_hashing=False)
    pairwise distance computation.
    """

    def __init__(self, use_hashing=False):
        self.use_hashing = use_hashing

    def fit(self, cdr3):
        if self.use_hashing:
            return to_dataframe(parallelized_clustering_on_length(cdr3, cluster_using_hash))
        else:
            return to_dataframe(parallelized_clustering_on_length(cdr3, cluster_using_distance))


def build_hash(cdr3):
    cdr3hash = {}
    for cdr in set(cdr3):
        for hash in (cdr[::2], cdr[1::2]):
            if hash not in cdr3hash:
                cdr3hash[hash] = set()
            cdr3hash[hash].add(cdr)
    return cdr3hash


def hd1(seq1, seq2):
    return sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2)) == 1


def edges_using_hash(cdr3):
    cdr3hash = build_hash(cdr3)
    edges = set()
    for hash in cdr3hash:
        if len(cdr3hash[hash]) < 2:
            continue
        for cdr1 in cdr3hash[hash]:
            for cdr2 in cdr3hash[hash]:
                if cdr1 != cdr2 and cdr1 <= cdr2 and hd1(cdr1, cdr2):
                    edges.add((cdr1, cdr2))
    return edges


def edges_using_distance(cdr3):
    edges = []
    for i in range(len(cdr3)):
        for j in range(i + 1, len(cdr3)):
            if hd1(cdr3[i], cdr3[j]):
                edges.append((cdr3[i], cdr3[j]))
    return edges


def cluster_using_distance(cdr3):
    return cluster(edges_using_distance(cdr3))


def cluster_using_hash(cdr3):
    return cluster(edges_using_hash(cdr3))


def parallelized_clustering_on_length(cdr3, clustering_func):
    grouped_on_length = defaultdict(list)
    for seq in cdr3:
        grouped_on_length[len(seq)].append(seq)
    grouped_on_length = grouped_on_length.values()
    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        res = parmap.map(clustering_func,
                         grouped_on_length,
                         pm_parallel=True,
                         pm_pool=pool)
    result = {}
    for i, clustering in enumerate(res):
        # Fix cluster ids
        for key in clustering:
            clustering[key] += i * 1000
        # Merge
        result.update(clustering)
    return result


def cluster(edges):
    clustering = {}
    cluster_id = 1
    for seq1, seq2 in edges:
        if seq1 not in clustering and seq2 not in clustering:
            clustering[seq1] = cluster_id
            clustering[seq2] = cluster_id
            cluster_id += 1
        elif seq1 not in clustering:
            clustering[seq1] = clustering[seq2]
        elif seq2 not in clustering:
            clustering[seq2] = clustering[seq1]
        elif clustering[seq1] != clustering[seq2]:
            cluster1 = clustering[seq1]
            cluster2 = clustering[seq2]
            for sequence in clustering:
                if clustering[sequence] == cluster2:
                    clustering[sequence] = cluster1
    return clustering


def to_dataframe(clustering):
    return pd.DataFrame({'CDR3': clustering.keys(), 'cluster': clustering.values()})


class Timer:
    def __init__(self, msg):
        self.msg = msg

    def __enter__(self):
        self.start = time.time()

    def __exit__(self, type, value, traceback):
        print(self.msg, time.time() - self.start)


def main():
    from clustcr import datasets, Clustering
    from clustcr.clustering.metrics import Metrics

    cdr3, epi = datasets.vdjdb_cdr3(), datasets.vdjdb_epitopes()

    with Timer('greedy'):
        clustering = GreedyClustering(use_hashing=True)
        result = clustering.fit(cdr3)
        print(result)
        print(Metrics(result, epi).summary())

    with Timer('normal'):
        print(Clustering(n_cpus='all').fit(cdr3).clusters_df)


if __name__ == '__main__':
    main()
