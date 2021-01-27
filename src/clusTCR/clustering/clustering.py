import pandas as pd
import multiprocessing


from .mcl import MCL, MCL_from_preclusters, MCL_multiprocessing_from_preclusters
from src.clusTCR.modules.faiss_clustering import FaissClustering
from .metrics import Metrics


class ClusteringResult:
    def __init__(self, nodelist):
        self.clusters_df = nodelist

    def metrics(self, epi):
        return Metrics(self.clusters_df, epi)


class Clustering:
    """
    The Clustering class offers flexible functionality for clustering of
    (large) sets of CDR3 amino acid sequences. The default clustering method
    of clusTCR is a two-step procedure that applies the faiss library to
    prune the search space and create superclusters. MCL is subsequently
    applied to identify specificity groups in each supercluster. This two-step
    clustering combination results in a fast and accurate way of grouping large
    data sets of CDR3 sequences into specificity groups.
    
    The Clustering module currently provides the following clustering methods:
        - MCL: Markov Clustering Algorithm. Accurate clustering method, 
        which is recommended for data sets containing < 50,000 sequences.
        - FAISS: This method provides extremely rapid clustering of sequences
        through dense vector representations. This method is far less accurate.
        This method also provides GPU support.
        - TWOSTEP: Combines the first two methods for a fast and accurate
        clustering method. This method provides both CPU multiprocessing,
        as well as GPU support.
    """

    def __init__(self,
                 method='two-step',
                 n_cpus=-1,
                 use_gpu=False,
                 faiss_cluster_size=5000,
                 mcl_params=None):

        """
        Parameters
        ----------
        cdr3 : pd.Series
            Input data consisting of a list of CDR3 amino acid sequences.
        method : str
            Clustering method. The default is two-step.
        n_cpus : int, optional
            Number of GPUs to use for clustering. -1 to use all GPUs.
            The default is -1.
        use_gpu : bool, optional
            Enable GPU computing for FAISS clustering. The default is False.
        mcl_params : list, optional
            Hyperparameters of MCL. The default is [1.2,2].
        faiss_cluster_size : TYPE, optional
            DESCRIPTION. The default is 5000.
        """
        self.mcl_params = mcl_params if mcl_params is not None else [1.2, 2]
        self.method = method.upper()
        self.use_gpu = use_gpu
        self.faiss_cluster_size = faiss_cluster_size
        self._set_n_cpus(n_cpus)

        available = ["MCL",
                     "FAISS",
                     "TWO-STEP"]
        assert self.method in available, f"Method not available, please choose one of the following methods:\n {available}"

    def _set_n_cpus(self, n_cpus):
        # Multiprocessing currently only available for Two-step method.
        # Use 1 cpu for other methods.
        if self.method != 'TWO-STEP' \
                or n_cpus < -1 \
                or n_cpus == 0 or \
                n_cpus > multiprocessing.cpu_count():
            self.n_cpus = 1
        elif n_cpus == -1:
            self.n_cpus = multiprocessing.cpu_count()
        else:
            self.n_cpus = n_cpus

    def _faiss(self, cdr3: pd.Series):
        """
        FAISS clustering method

        Returns
        -------
        clusters : pd.DataFrame
            pd.DataFrame containing two columns: 'CDR3' and 'cluster'.
            The first column contains CDR3 sequences, the second column
            contains the corresponding cluster ids.
        """
        clusters = {"CDR3": [], "cluster": []}
        faiss_output = FaissClustering.cluster(cdr3.reset_index(drop=True),
                                               avg_items_per_cluster=self.faiss_cluster_size,
                                               use_gpu=self.use_gpu)
        for cluster in faiss_output.get_cluster_contents():
            clusters["CDR3"].append(cluster)
            clusters["CDR3"].append(len(clusters["cluster"]))

        return pd.DataFrame(clusters)

    def _twostep(self, cdr3):
        """
        Two-step clustering procedure for speeding up CDR3 clustering by
        pre-sorting sequences into superclusters. A second clustering step
        is performed on each individual supercluster.

        Parameters
        ----------
        supercluster_size : int, optional
            Approximate number of sequences in each supercluster. 
            The default is 5000.

        Returns
        -------
        nodelist : pd.DataFrame
            pd.DataFrame containing two columns: 'CDR3' and 'cluster'.
            The first column contains CDR3 sequences, the second column
            contains the corresponding cluster ids.
        """

        super_clusters = FaissClustering.cluster(cdr3.reset_index(drop=True),
                                                 avg_items_per_cluster=self.faiss_cluster_size,
                                                 use_gpu=self.use_gpu)
        if self.n_cpus > 1:
            return MCL_multiprocessing_from_preclusters(cdr3, super_clusters, self.n_cpus)
        else:
            return MCL_from_preclusters(cdr3, super_clusters)

    def fit(self, cdr3: pd.Series):
        """
        Function that calls the indicated clustering method and returns clusters
        in a nodelist format.
        
        Parameters
        ----------

        Returns
        -------
        nodelist : TYPE
            Table containing sequences and their corresponding cluster ids.
        """
        if self.method == 'MCL':
            return ClusteringResult(MCL(cdr3, mcl_hyper=self.mcl_params))
        elif self.method == 'FAISS':
            return ClusteringResult(self._faiss(cdr3))
        else:
            return ClusteringResult(self._twostep(cdr3))
