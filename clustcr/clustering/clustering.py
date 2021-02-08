import pandas as pd
import multiprocessing
from typing import Union
from os.path import join, exists
from os import mkdir, getcwd
from shutil import rmtree
import random

from .mcl import MCL, MCL_from_preclusters, MCL_multiprocessing_from_preclusters
from clustcr.modules.faiss_clustering import FaissClustering
from clustcr.analysis.features import FeatureGenerator
from .metrics import Metrics


class ClusteringResult:
    def __init__(self, nodelist):
        self.clusters_df = nodelist
        
    def summary(self):
        motifs = FeatureGenerator(self.clusters_df).clustermotif()
        summ = self.clusters_df.cluster.value_counts().to_frame()
        summ.rename(columns={'cluster':'size'},inplace=True)
        summ = summ.rename_axis('cluster_idx').reset_index()
        summ['motif'] = motifs.values()
        return summ
    
    def write_to_csv(self, path=join(getcwd(),'clusTCR_clusters.csv')):
        return self.clusters_df.to_csv(path,index=False)

    def cluster_contents(self):
        return list(self.clusters_df.groupby(['cluster'])['CDR3'].apply(list))

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

    BATCH_TMP_DIRECTORY = 'clustcr_batch_tmp' + str(random.randint(0, 10 ** 8))

    def __init__(self,
                 method='two-step',
                 n_cpus: Union[str, int] = 'all',
                 use_gpu=False,
                 faiss_cluster_size=5000,
                 mcl_params=None,
                 faiss_training_data=None,
                 max_sequence_size=None,
                 fitting_data_size=None):

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

        # For batch processing
        self.faiss_training_data = faiss_training_data
        self.max_sequence_size = max_sequence_size
        if fitting_data_size and self.faiss_training_data is not None:
            self.faiss_cluster_size = int(self.faiss_cluster_size / (fitting_data_size / len(self.faiss_training_data)))
            self.faiss_clustering = self._train_faiss(faiss_training_data)
            if exists(Clustering.BATCH_TMP_DIRECTORY):
                rmtree(Clustering.BATCH_TMP_DIRECTORY)
            mkdir(Clustering.BATCH_TMP_DIRECTORY)
        else:
            self.faiss_clustering = None

        self._set_n_cpus(n_cpus)
        available = ["MCL",
                     "FAISS",
                     "TWO-STEP"]
        assert self.method in available, f"Method not available, please choose one of the following methods:\n {available}"

    def _set_n_cpus(self, n_cpus):
        # Multiprocessing currently only available for Two-step method.
        # Use 1 cpu for other methods.
        if self.method != 'TWO-STEP' \
                or (isinstance(n_cpus, str) and n_cpus != 'all') \
                or (isinstance(n_cpus, int) and (n_cpus < 1 or n_cpus > multiprocessing.cpu_count())):
            self.n_cpus = 1
        elif n_cpus == 'all':
            self.n_cpus = multiprocessing.cpu_count()
        else:
            self.n_cpus = n_cpus

    def _train_faiss(self, cdr3: pd.Series):
        clustering = FaissClustering(avg_cluster_size=self.faiss_cluster_size,
                                     use_gpu=self.use_gpu,
                                     max_sequence_size=self.max_sequence_size)
        clustering.train(cdr3)
        return clustering

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
        cdr3 = cdr3.reset_index(drop=True)
        if self.faiss_clustering is not None:
            clustering = self.faiss_clustering
        else:
            clustering = self._train_faiss(cdr3)

        result = clustering.cluster(cdr3)
        clusters = {"CDR3": [], "cluster": []}
        for i, cluster in enumerate(result):
            clusters["CDR3"].append(cdr3[i])
            clusters["cluster"].append(int(cluster))

        return ClusteringResult(pd.DataFrame(clusters))

    def _twostep(self, cdr3):
        """
        Two-step clustering procedure for speeding up CDR3 clustering by
        pre-sorting sequences into superclusters. A second clustering step
        is performed on each individual supercluster.

        Parameters
        ----------
        cdr3 : pd.Series
            Input data consisting of a list of CDR3 amino acid sequences.

        Returns
        -------
        nodelist : pd.DataFrame
            pd.DataFrame containing two columns: 'CDR3' and 'cluster'.
            The first column contains CDR3 sequences, the second column
            contains the corresponding cluster ids.
        """

        super_clusters = self._faiss(cdr3)
        if self.n_cpus > 1:
            return ClusteringResult(MCL_multiprocessing_from_preclusters(cdr3, super_clusters, self.n_cpus))
        else:
            return ClusteringResult(MCL_from_preclusters(cdr3, super_clusters))

    def batch_precluster(self, cdr3: pd.Series):
        assert self.faiss_clustering is not None, 'Batch precluster needs faiss_training_data and fitting_data_size'
        clustered = self._faiss(cdr3)
        for index, row in clustered.clusters_df.iterrows():
            filename = join(Clustering.BATCH_TMP_DIRECTORY, str(row['cluster']))
            with open(filename, 'a') as f:
                f.write(row['CDR3'] + '\n')

    def batch_cluster(self):
        """
        Clusters the preclusters (stored on disk) using MCL.
        Thus requires the batch_precluster method to be called beforehand.

        - Multiprocessing is used (if enabled) to cluster multiple preclusters at the same time.
        - Generator function (by using yield) to prevent memory overflow

        The amount of preclusters that is clustered by MCL per batch (iteration) is calculated as such:
            - Check how many preclusters roughly contain 50k sequences when combined (as that should fit in memory no problem)
            - Limit to bounds (1, ncpus)
        """
        clusters_per_batch = max(1, min(self.n_cpus, 50000 // self.faiss_cluster_size))
        npreclusters = self.faiss_clustering.ncentroids()
        max_cluster_id = 0
        for i in range(0, npreclusters, clusters_per_batch):
            cluster_ids = range(i, min(i + clusters_per_batch, npreclusters))
            preclusters = self._batch_process_preclusters(cluster_ids)
            mcl_result = MCL_multiprocessing_from_preclusters(None, preclusters, self.n_cpus)
            mcl_result['cluster'] += max_cluster_id + 1
            max_cluster_id = mcl_result['cluster'].max()
            yield ClusteringResult(mcl_result)

    def _batch_process_preclusters(self, cluster_ids):
        preclusters = {'CDR3': [], 'cluster': []}
        for cluster_id in cluster_ids:
            filename = join(Clustering.BATCH_TMP_DIRECTORY, str(cluster_id))
            if not exists(filename):
                continue
            with open(filename) as f:
                sequences = f.readlines()
                preclusters['CDR3'].extend(sequences)
                preclusters['cluster'].extend([cluster_id] * len(sequences))
        return ClusteringResult(pd.DataFrame(preclusters))

    def batch_cleanup(self):
        rmtree(Clustering.BATCH_TMP_DIRECTORY)

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
            return self._faiss(cdr3)
        else:
            return self._twostep(cdr3)
