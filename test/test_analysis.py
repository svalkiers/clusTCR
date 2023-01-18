from base import TestBase
from clustcr import datasets, Clustering, ClusterAnalysis, ModelTraining


class ClusteringTest(TestBase):

    def setUp(self):
        self.cdr3 = datasets.test_cdr3().junction_aa
        self.epitopes = datasets.test_epitopes()[["junction_aa", "epitope"]]
        self.clustering_result = Clustering().fit(self.cdr3)

    def make_features(self):
        return self.clustering_result.compute_features(compute_pgen=True)

    def test_feature_generation(self):
        self.make_features()

    def test_pca(self):
        ClusterAnalysis(self.make_features()).pca()

    def test_prediction(self):
        ClusterAnalysis(self.make_features()).predict_quality()

#    def test_train_model(self):
#        model = ModelTraining(self.clustering_result.clusters_df, self.epitopes)
#        fitted = model.fit_data()
#        model.evaluate()
#        model.save(fitted, 'test.pkl')


