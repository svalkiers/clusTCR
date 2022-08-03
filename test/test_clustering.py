from base import TestBase
from clustcr import datasets, Clustering


class ClusteringTest(TestBase):

    def setUp(self):
        self.cdr3 = datasets.test_cdr3()["junction_aa"]
        self.epitopes = datasets.test_epitopes()[["junction_aa", "epitope"]]

    def test_normal(self):
        Clustering().fit(self.cdr3)

    def test_quality(self):
        metrics = Clustering().fit(self.cdr3).metrics(self.epitopes)
        self.assertGreater(metrics.purity()[0], 0.6)
        self.assertGreater(metrics.consistency()[0], 0.12)
        self.assertGreater(metrics.retention(), 0.21)
        self.assertGreater(metrics.purity_90()[0], 0.36)

    def test_mcl(self):
        Clustering(method='mcl').fit(self.cdr3)

    def test_faiss(self):
        Clustering(method='faiss').fit(self.cdr3)

    def test_alphabeta(self):
        df = datasets.vdjdb_paired()
        alpha, beta = df['CDR3_alpha'], df['CDR3_beta']
        Clustering().fit(beta, alpha=alpha)

    def test_multiprocessing(self):
        for cpu in [-1, 0, 1, 2, 'all']:
            for method in ['two-step', 'faiss', 'mcl']:
                Clustering(method=method, n_cpus=cpu).fit(self.cdr3)

    def test_faiss_cluster_size(self):
        for size in range(2, 6003, 2000):
            for method in ['two-step', 'faiss', 'mcl']:
                Clustering(method=method, faiss_cluster_size=size).fit(self.cdr3)

    def test_batch_clustering(self):
        vdj = datasets.vdjdb_beta()
        max_sequence_size = vdj.junction_aa.str.len().max()
        train = vdj.sample(2000)
        times = 3
        size_per_time = 3000
        clustering = Clustering(faiss_training_data=train.junction_aa,
                                fitting_data_size=times * size_per_time,
                                max_sequence_size=max_sequence_size)
        for i in range(times):
            sample = vdj.sample(size_per_time)
            clustering.batch_precluster(sample.junction_aa)
        for clusters in clustering.batch_cluster():
            df = clusters.clusters_df
        clustering.batch_cleanup()

    def test_batch_clustering_multiprocessing(self):
        vdj = datasets.vdjdb_beta()
        max_sequence_size = vdj.junction_aa.str.len().max()
        train = vdj.sample(2000)
        times = 3
        size_per_time = 3000
        clustering = Clustering(faiss_training_data=train.junction_aa,
                                fitting_data_size=times * size_per_time,
                                max_sequence_size=max_sequence_size,
                                n_cpus='all')
        for i in range(times):
            sample = vdj.sample(size_per_time)
            clustering.batch_precluster(sample.junction_aa)
        for clusters in clustering.batch_cluster():
            df = clusters.clusters_df
        clustering.batch_cleanup()

    def test_matrix(self):
        vdj = datasets.vdjdb_beta()
        max_sequence_size = vdj.junction_aa.str.len().max()
        train = vdj.sample(2000)
        times = 3
        size_per_time = 3000
        clustering = Clustering(faiss_training_data=train.junction_aa,
                                fitting_data_size=times * size_per_time,
                                max_sequence_size=max_sequence_size)
        for i in range(times):
            sample = vdj.sample(size_per_time)
            clustering.batch_precluster(sample.junction_aa, name=f'time {i}')
        for clusters in clustering.batch_cluster(calc_cluster_matrix=True):
            df = clusters.clusters_df
        clustering.batch_cluster_matrix()
        clustering.batch_cleanup()

    def test_metrics(self):
        metrics = Clustering().fit(self.cdr3).metrics(self.epitopes)
        metrics.purity()
        metrics.consistency()
        metrics.retention()
        metrics.purity_90()
        metrics.summary()

    def test_summary(self):
        Clustering().fit(self.cdr3).summary()

    def test_write_to_csv(self):
        Clustering().fit(self.cdr3).write_to_csv()

    def test_cluster_contents(self):
        Clustering().fit(self.cdr3).cluster_contents()


