from base import TestBase
from clustcr import datasets, Clustering


class ClusteringTest(TestBase):

    def setUp(self):
        self.cdr3 = datasets.test_cdr3()
        self.epitopes = datasets.test_epitopes()

    def test_normal(self):
        Clustering().fit(self.cdr3)

    def test_mcl(self):
        Clustering(method='mcl').fit(self.cdr3)

    def test_faiss(self):
        Clustering(method='faiss').fit(self.cdr3)

    def test_multiprocessing(self):
        for method in ['two-step', 'faiss', 'mcl']:
            Clustering(method=method, n_cpus='all').fit(self.cdr3)

    def test_faiss_cluster_size(self):
        for size in range(2, 6003, 2000):
            for method in ['two-step', 'faiss', 'mcl']:
                Clustering(method=method, faiss_cluster_size=size).fit(self.cdr3)

    def test_metrics(self):
        metrics = Clustering().fit(self.cdr3).metrics(self.epitopes)
        metrics.purity()
        metrics.consistency()
        metrics.retention()
        metrics.purity_90()

