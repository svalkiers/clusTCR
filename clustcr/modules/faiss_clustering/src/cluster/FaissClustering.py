import pandas as pd
import faiss

from ..profile.properties import *
from .tools import make_profiles


class FaissClustering:

    def __init__(self, max_sequence_size=None, properties: list = OPTIMAL, avg_cluster_size=10, use_gpu=False):
        self.properties = properties
        self.avg_cluster_size = avg_cluster_size
        self.use_gpu = use_gpu
        self.kmeans = None
        self.max_sequence_size = max_sequence_size

    def train(self,
              data: pd.Series,
              properties: list = OPTIMAL,
              avg_cluster_size=10,
              use_gpu=False):
        """
        Clusters the given CDR3 series and returns a FaissClustering object

        @params
        properties: to be used during the clustering can be given in a list, and are available to import as
                        'from faiss_clustering.properties import *'

        avg_cluster_size: the average size of the clusters that will be generated
        """
        profiles = make_profiles(data, properties, self.max_sequence_size)
        size, dimension = profiles.shape
        centroids = max(1, size // avg_cluster_size)
        self.kmeans = faiss.Kmeans(dimension, centroids, gpu=use_gpu, min_points_per_centroid=1)
        self.kmeans.train(profiles)

    def cluster(self, data: pd.Series):
        profiles = make_profiles(data, self.properties, self.max_sequence_size)
        D, I = self.kmeans.index.search(profiles, 1)
        return I
