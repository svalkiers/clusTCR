import multiprocessing
import parmap
import pandas as pd
from collections import defaultdict
import itertools


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
        self.cdr3hash = None

    def fit(self, cdr3):
        if self.use_hashing:
            return GreedyClustering.to_dataframe(self.cluster(self.edges_from_hash(cdr3)))
        else:
            return GreedyClustering.to_dataframe(self.parallelized_clustering_using_distance(cdr3))

    def build_hash(self, cdr3):
        self.cdr3hash = {}
        for cdr in set(cdr3):
            for hash in (cdr[::2], cdr[1::2]):
                if hash not in self.cdr3hash:
                    self.cdr3hash[hash] = set()
                self.cdr3hash[hash].add(cdr)

    @staticmethod
    def hd1(seq1, seq2):
        return sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2)) == 1

    def edges_from_hash(self, cdr3):
        self.build_hash(cdr3)
        edges = set()
        for hash in self.cdr3hash:
            if len(self.cdr3hash[hash]) < 2:
                continue
            for cdr1 in self.cdr3hash[hash]:
                for cdr2 in self.cdr3hash[hash]:
                    if cdr1 != cdr2 and cdr1 <= cdr2 and GreedyClustering.hd1(cdr1, cdr2):
                        edges.add((cdr1, cdr2))
        return edges

    @staticmethod
    def edges_using_distance(cdr3):
        edges = []
        for i in range(len(cdr3)):
            for j in range(i + 1, len(cdr3)):
                if GreedyClustering.hd1(cdr3[i], cdr3[j]):
                    edges.append((cdr3[i], cdr3[j]))
        return edges

    @staticmethod
    def cluster_using_distance(cdr3):
        return GreedyClustering.cluster(GreedyClustering.edges_using_distance(cdr3))

    @staticmethod
    def parallelized_clustering_using_distance(cdr3):
        grouped_on_length = defaultdict(list)
        for seq in cdr3:
            grouped_on_length[len(seq)].append(seq)
        grouped_on_length = grouped_on_length.values()
        with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
            res = parmap.map(GreedyClustering.cluster_using_distance,
                             grouped_on_length,
                             pm_parallel=True,
                             pm_pool=pool)
        result = {}
        for i, clustering in enumerate(res):
            for key in clustering:
                clustering[key] += i * 1000
            result.update(clustering)
        return result

    @staticmethod
    def edges_from_parallelization(cdr3):
        grouped_on_length = defaultdict(list)
        for seq in cdr3:
            grouped_on_length[len(seq)].append(seq)
        grouped_on_length = grouped_on_length.values()
        with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
            res = parmap.map(GreedyClustering.edges_using_distance,
                             grouped_on_length,
                             pm_parallel=True,
                             pm_pool=pool)
        return itertools.chain.from_iterable(res)

    @staticmethod
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

    @staticmethod
    def to_dataframe(clustering):
        return pd.DataFrame({'CDR3': clustering.keys(), 'cluster': clustering.values()})


if __name__ == '__main__':
    from clustcr import datasets
    from clustcr.clustering.metrics import Metrics

    cdr3, epi = datasets.vdjdb_cdr3(), datasets.vdjdb_epitopes()
    clustering = GreedyClustering(use_hashing=True)
    result = clustering.fit(cdr3)
    print(result)
    print(Metrics(result, epi).summary())
