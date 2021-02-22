from base import TestBase
from clustcr import datasets, Clustering


class ClusteringTest(TestBase):

    def setUp(self):
        self.cdr3 = datasets.test_cdr3()
        self.epitopes = datasets.test_epitopes()

    def test_normal(self):
        Clustering().fit(self.cdr3)

    def test_metrics(self):
        metrics = Clustering().fit(self.cdr3).metrics()
        metrics.purity()
        metrics.consistency()
        metrics.retention()
        metrics.purity_90()

