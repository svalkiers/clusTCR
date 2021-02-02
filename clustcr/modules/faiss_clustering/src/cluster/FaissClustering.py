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

    def train(self, data: pd.Series):
        """
        Trains the kmeans clustering using the given cdr3 data
        """
        profiles = make_profiles(data, self.properties, self.max_sequence_size)
        size, dimension = profiles.shape
        centroids = max(1, size // self.avg_cluster_size)
        self.kmeans = faiss.Kmeans(dimension, centroids, gpu=self.use_gpu, min_points_per_centroid=1)
        self.kmeans.train(profiles)

    def cluster(self, data: pd.Series):
        assert self.kmeans is not None, "FaissClustering not trained"
        profiles = make_profiles(data, self.properties, self.max_sequence_size)
        D, I = self.kmeans.index.search(profiles, 1)
        return I
