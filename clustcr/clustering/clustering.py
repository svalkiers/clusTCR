import pandas as pd
import multiprocessing
from typing import Union
from os.path import join, exists
from os import mkdir, getcwd
from shutil import rmtree
import random
import time

from .mcl import MCL, MCL_from_preclusters, MCL_multiprocessing_from_preclusters
from clustcr.modules.faiss_clustering import FaissClustering, properties
from clustcr.analysis.features import FeatureGenerator
from .metrics import Metrics
from .multirepertoire_cluster_matrix import MultiRepertoireClusterMatrix
from .tools import create_edgelist

def timeit(myfunc):
    # Decorator to keep track of time required to run a function
    def timed(*args, **kwargs):
        start = time.time()
        result = myfunc(*args, **kwargs)
        end = time.time()
        print(f'Total time to run \'{myfunc.__name__}\': {(end-start):.3f}s')
        return result
    return timed

class ClusteringResult:
    def __init__(self, nodelist):
        self.clusters_df = nodelist

    def summary(self, motif_cutoff=.7):
        motifs = FeatureGenerator(self.clusters_df).clustermotif(cutoff=motif_cutoff)
        summ = self.clusters_df.cluster.value_counts().to_frame()
        summ.rename(columns={'cluster': 'size'}, inplace=True)
        summ = summ.rename_axis('cluster_idx').reset_index()
        summ['motif'] = motifs.values()
        return summ

    def write_to_csv(self, path=join(getcwd(), 'clusTCR_clusters.csv')):
        return self.clusters_df.to_csv(path, index=False)

    def export_network(self, filename='clusTCR_network.txt'):
        return create_edgelist(self.clusters_df.CDR3, filename)

    def cluster_contents(self):
        return list(self.clusters_df.groupby(['cluster'])['CDR3'].apply(list))

    def compute_features(self, compute_pgen=True):
        return FeatureGenerator(self.clusters_df).get_features(compute_pgen=compute_pgen)

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
                 distance_metric='hamming',
                 n_cpus: Union[str, int] = 1,
                 use_gpu=False,
                 faiss_cluster_size=5000,
                 mcl_params=None,
                 faiss_training_data=None,
                 max_sequence_size=None,
                 fitting_data_size=None,
                 rnd_chunk_size=5000):

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
        self.distance_metric = distance_metric.upper()
        self.use_gpu = use_gpu
        self.faiss_cluster_size = faiss_cluster_size
        self.faiss_properties = properties.OPTIMAL
        self._set_n_cpus(n_cpus)
        self.rnd_chunk_size = rnd_chunk_size

        # For batch processing
        self.faiss_training_data = faiss_training_data
        self.max_sequence_size = max_sequence_size
        self.cluster_matrix = None
        if fitting_data_size and self.faiss_training_data is not None:
            self.faiss_cluster_size = int(self.faiss_cluster_size / (fitting_data_size / len(self.faiss_training_data)))
            self.faiss_clustering = self._train_faiss(faiss_training_data)
            if exists(Clustering.BATCH_TMP_DIRECTORY):
                rmtree(Clustering.BATCH_TMP_DIRECTORY)
            mkdir(Clustering.BATCH_TMP_DIRECTORY)
        else:
            self.faiss_clustering = None

        available = ["MCL",
                     "FAISS",
                     "TWO-STEP",
                     "RANDOM"]
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

    def _random(self, cdr3):
        n = round(len(cdr3) / self.rnd_chunk_size)
        sequences = list(cdr3)
        random.shuffle(sequences)
        chunks = [sequences[i::n] for i in range(n)]
        clusters = {"CDR3": [], "cluster": []}
        for i, chunk in enumerate(chunks):
            for seq in chunk:
                clusters["CDR3"].append(seq)
                clusters["cluster"].append(i)
        return ClusteringResult(pd.DataFrame(clusters))

    def _train_faiss(self, cdr3: pd.Series, get_profiles=False):
        clustering = FaissClustering(avg_cluster_size=self.faiss_cluster_size,
                                     use_gpu=self.use_gpu,
                                     max_sequence_size=self.max_sequence_size,
                                     properties=self.faiss_properties,
                                     n_cpus=self.n_cpus)
        profiles = clustering.train(cdr3)
        if get_profiles:
            return clustering, profiles
        else:
            return clustering

    def _faiss(self, cdr3: pd.Series) -> ClusteringResult:
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
        profiles = None
        if self.faiss_clustering is not None:
            clustering = self.faiss_clustering
        else:
            clustering, profiles = self._train_faiss(cdr3, get_profiles=True)

        if profiles is not None:
            result = clustering.cluster(profiles, is_profile=True)
        else:
            result = clustering.cluster(cdr3)

        clusters = {"CDR3": [], "cluster": []}
        for i, cluster in enumerate(result):
            clusters["CDR3"].append(cdr3[i])
            clusters["cluster"].append(int(cluster))

        return ClusteringResult(pd.DataFrame(clusters))

    @timeit
    def _twostep(self, cdr3) -> ClusteringResult:
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
            return ClusteringResult(MCL_multiprocessing_from_preclusters(cdr3, super_clusters, self.distance_metric, self.mcl_params, self.n_cpus))
        else:
            return ClusteringResult(MCL_from_preclusters(cdr3, super_clusters, self.distance_metric, self.mcl_params))

    def batch_precluster(self, cdr3: pd.Series, name=''):
        assert self.faiss_clustering is not None, 'Batch precluster needs faiss_training_data and fitting_data_size'
        clustered = self._faiss(cdr3)
        for index, row in clustered.clusters_df.iterrows():
            filename = join(Clustering.BATCH_TMP_DIRECTORY, str(row['cluster']))
            with open(filename, 'a') as f:
                f.write(f'{row["CDR3"]},{name}\n')

    def batch_cluster(self, calc_cluster_matrix=False):
        """
        Clusters the preclusters (stored on disk) using MCL.
        Thus requires the batch_precluster method to be called beforehand.

        - Multiprocessing is used (if enabled) to cluster multiple preclusters at the same time.
        - Generator function (by using yield) to prevent memory overflow

        The amount of preclusters that is clustered by MCL per batch (iteration) is calculated as such:
            - Check how many preclusters roughly contain 50k sequences when combined (as that should fit in memory no problem)
            - Limit to bounds (1, ncpus)
        """
        self.cluster_matrix = MultiRepertoireClusterMatrix()
        clusters_per_batch = max(1, min(self.n_cpus, 50000 // self.faiss_cluster_size))
        npreclusters = self.faiss_clustering.ncentroids()
        max_cluster_id = 0
        for i in range(0, npreclusters, clusters_per_batch):
            cluster_ids = range(i, min(i + clusters_per_batch, npreclusters))
            preclusters = self._batch_process_preclusters(cluster_ids)
            mcl_result = MCL_multiprocessing_from_preclusters(None, preclusters, self.distance_metric, self.mcl_params, self.n_cpus)
            mcl_result['cluster'] += max_cluster_id + 1
            max_cluster_id = mcl_result['cluster'].max()
            if calc_cluster_matrix:
                self.cluster_matrix.add(preclusters, mcl_result)
            yield ClusteringResult(mcl_result)

    def _batch_process_preclusters(self, cluster_ids):
        preclusters = {'CDR3': [], 'cluster': [], 'name': []}
        for cluster_id in cluster_ids:
            filename = join(Clustering.BATCH_TMP_DIRECTORY, str(cluster_id))
            if not exists(filename):
                continue
            with open(filename) as f:
                content = f.readlines()
                content = [e.replace('\n', '').split(',') for e in content]
                sequences = list(map(lambda x: x[0], content))
                names = list(map(lambda x: x[1], content))
                preclusters['CDR3'].extend(sequences)
                preclusters['cluster'].extend([cluster_id] * len(sequences))
                preclusters['name'].extend(names)
        return ClusteringResult(pd.DataFrame(preclusters))

    def batch_cluster_matrix(self):
        return self.cluster_matrix.get_matrix()

    def batch_cleanup(self):
        rmtree(Clustering.BATCH_TMP_DIRECTORY)

    def fit(self, cdr3: pd.Series, alpha: pd.Series = None) -> ClusteringResult:
        """
        Function that calls the indicated clustering method and returns clusters in a ClusteringResult
        """
        if alpha is not None:
            assert len(cdr3) == len(alpha), 'amount of CDR3 data is not equal to amount of alpha chain data'
            cdr3 = cdr3.add(alpha)
        if self.method == 'MCL':
            return ClusteringResult(MCL(cdr3, distance_metric=self.distance_metric, mcl_hyper=self.mcl_params))
        elif self.method == 'FAISS':
            return self._faiss(cdr3)
        elif self.method == 'RANDOM':
            return self._random(cdr3)
        else:
            return self._twostep(cdr3)
