import numpy as np
import pandas as pd
import faiss

from ..profile.properties import *
from .tools import make_profiles, cluster_with_faiss


class FaissClustering:

    @staticmethod
    def cluster(data: pd.Series,
                properties: list = OPTIMAL,
                avg_items_per_cluster=10,
                use_gpu=False,
                vector_mapping_func=None,
                add_average=False,
                add_length=False):
        """
        Clusters the given CDR3 series and returns a FaissClustering object

        @params
        properties: to be used during the clustering can be given in a list, and are available to import as
                        'from faiss_clustering.properties import *'

        avg_items_per_cluster: the average size of the clusters that will be generated

        vector_mapping_func: function that takes a profile (a list of float values representing the sequence)
                                and returns a list with the same size

        add_average: adds the average of each of the wanted properties over a sequences to the profile vector
        add_length: adds the length of the sequence to the profile
        """
        profiles = make_profiles(data, properties,
                                 vector_mapping_func=vector_mapping_func,
                                 add_average=add_average,
                                 add_length=add_length)
        clustering = cluster_with_faiss(profiles,
                                        items_per_cluster=avg_items_per_cluster,
                                        ids=data.keys().to_numpy(),
                                        use_gpu=use_gpu)
        if use_gpu:
            clustering = clustering.index_
        return FaissClustering(clustering, data)

    def __init__(self, index, data):
        self.index = index
        self.data = data

    def get_cluster_contents(self, include_sequences=True):
        """
        Returns a list that contains all the clusters
        Each cluster is represented as a list
        For example:
            [[CAI, CAT, CTT], [TAA, CAA]]
        where there are two clusters: (CAI, CAT, CTT) and (TAA, CAA)

        include_sequences:
            if True, will include the full sequences as shown in the examples
            if False, will include IDs of the sequences
        """
        cluster_contents = []
        for i in range(self.index.nlist):
            size = self.index.invlists.list_size(i)
            ids = list(faiss.rev_swig_ptr(self.index.invlists.get_ids(i), size))
            if include_sequences:
                cluster_contents.append([self.data[id] for id in ids])
            else:
                cluster_contents.append(ids)
        return cluster_contents
