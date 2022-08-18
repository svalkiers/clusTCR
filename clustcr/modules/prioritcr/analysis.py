from collections.abc import Sequence
from collections import Counter
from functools import cached_property, partial
import io
from itertools import count, chain
from typing import List, Union
from multiprocessing import Pool

from clustcr.clustering.clustering import ClusteringResult
from datasketch import MinHash, MinHashLSH, MinHashLSHEnsemble, MinHashLSHForest
from logomaker import Logo, transform_matrix
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .generation_probability import PGen
from .physicochemical_properties import AALPHABET, AA_DICT, BG_POSITION_MATRICES
from .tools import kmer_iterator, hash_kmer, vec_bin_array


class Cluster(Sequence):
    id_iter = count()  # counter for ids
    """
    Class for storing a cluster (=list of CDR3 sequences) and calculate
    its useful properties.
    """

    def __init__(
        self,
        sequences: List[str],
        cluster_id: int = None,
        k: Union[int, range, list, str] = 5,
    ) -> None:
        self.sequences = sequences
        self.xid = cluster_id if cluster_id else next(Cluster.id_iter)
        self._k = k
        self._k_list = self.kmerize(k)
        self.cluster_size = len(self.sequences)
        self.regex = self.create_regex()
        self._num_perm = 128

    def __repr__(self) -> str:
        return f"id={self.xid}, cluster_size={self.cluster_size}, sequence_len={self.sequence_len}, regex={self.regex}"

    def __len__(self) -> int:
        return self.cluster_size

    def __getitem__(self, i):
        return self.sequences[i]

    def __iter__(self):
        return (x for x in self.sequences)

    def __add__(self, other):
        assert self.k == other.k
        return Cluster(list(set(chain(self, other))), k=self.k)

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return other.__add__(self)

    def __eq__(self, other) -> bool:
        return self.xid == other.xid

    def __hash__(self) -> int:
        return hash(f"{self.xid}{self.regex}")

    def _create_profile_matrix(self, return_df: bool = True):
        # only proceed if all sequences have the same lenght, else remove seq
        seq_list = [s for s in self.sequences if len(s) == self.sequence_len]

        # initiate profile matrix
        pm = np.zeros(shape=(len(AALPHABET), self.sequence_len))

        # initiate AA dict:
        AAs = {AA: i for i, AA in enumerate(AALPHABET)}

        # update profile matrix with counts
        for s in seq_list:
            for i, AA in enumerate(s):
                pm[AAs[AA], i] += 1

        # normalize profile matrix to percentages
        pm = pm / len(seq_list)

        if return_df:
            # convert numpy result matrix to dataframe
            cols = [f"p{i}" for i in range(self.sequence_len)]
            pm = pd.DataFrame(pm, columns=cols, index=list(AALPHABET))

        return pm

    @property
    def sequence_len(self):
        """
        Most common length of CDR3 sequence in cluster.
        """
        lenlist = [len(x) for x in self.sequences]
        return Counter(lenlist).most_common()[0][0]

    @property
    def k(self):
        """
        Value used to split cluster sequences into k-mers.

        Additionally; a list or range of values, or 'all' can be passed to use multiple k values.
        """
        return self._k

    @k.setter
    def k(self, k: Union[int, range, list, str]):
        if k != self._k:
            self._k = k
            self._k_list = self.kmerize(k)
            self._reset_cached_properties()

    def _reset_cached_properties(self):
        for p in ["mhash", "shash_32", "shash_64", "shash_128"]:
            if p in self.__dict__:
                del self.__dict__[p]

    @cached_property
    def mhash(self):
        """
        MinHash hash of cluster kmers.
        """
        m = MinHash(self.num_perm)
        for s in self._k_list:
            m.update(s.encode("utf8"))
        return m

    @cached_property
    def shash_32(self):
        """
        32-bit simhash of cluster kmers.
        """
        return self._shash(m=32)

    @cached_property
    def shash_64(self):
        """
        64-bit simhash of cluster kmers.
        """
        return self._shash(m=64)

    @cached_property
    def shash_128(self):
        """
        128-bit simhash of cluster kmers.
        """
        return self._shash(m=128)

    def _shash(self, m):
        """
        Simhash creation function
        """
        total_kmer_count = len(self._k_list)
        weight_dict = {
            k: v / total_kmer_count for k, v in Counter(self._k_list).items()
        }
        kmers = list(set(self._k_list))
        weights = np.array([weight_dict[km] for km in kmers])
        hashed_kmers = [hash_kmer(k, m) for k in kmers]
        binr = np.array(vec_bin_array(hashed_kmers, m), dtype=np.float16)
        binr[binr == 0] = -1
        w = (binr.T * weights).T
        w_s = w.sum(axis=0)
        simh = np.where(w_s < 0, 0, 1)
        return simh

    @property
    def num_perm(self):
        return self._num_perm

    @num_perm.setter
    def num_perm(self, num_perm):
        if "mhash" in self.__dict__:
            del self.__dict__["mhash"]
        self._num_perm = num_perm

    @cached_property
    def pgens(self):
        """
        A pd.Series containing the generation probabilities of cluster sequences.
        """
        pgens = list(PGen(parallel=False).compute_multiple(self.sequences))
        return pd.Series(pgens, index=self.sequences)

    @property
    def average_pgen(self) -> float:
        """
        Average generation probability of cluster sequences.
        """
        return self.pgens.mean()

    @property
    def average_surprise(self) -> float:
        """
        Average -10log of the generation probability of cluster sequences.
        """
        s = -np.log10(self.pgens)
        return s.mean()

    @property
    def read_count(self) -> int:
        """
        Sum of read counts of all CDR3 sequences in the cluster

        read_count can be set from an int but also from a counter or dict with
        the CDR3 sequence as keys, and the read count as values.
        """
        return self._read_count

    @read_count.setter
    def read_count(self, reads: Union[int, Counter, dict]) -> None:
        if isinstance(reads, int):
            self._read_count = reads
        else:
            self._read_count = sum([reads.get(s, 0) for s in self.sequences])

    def profile_matrix(self):
        return self._create_profile_matrix()

    def create_regex(self, min_cutoff: float = 0.1) -> str:
        """
        Create regex representation of the cluster.

        Parameters
        ----------
        cutoff : float, default=0.1
            The minimum probability of an AA at a position to be included in the regex.
        """
        pm = self._create_profile_matrix(return_df=False)
        res = ""
        for pos in range(self.sequence_len):
            aa_options = np.where(pm[:, pos] > min_cutoff)[0]
            if len(aa_options) == 1:
                res += AA_DICT[aa_options[0]]
            elif len(aa_options) == 0:
                res += "."
            else:
                res += f'[{"".join([AA_DICT[x] for x in aa_options])}]'
        return res

    def plot_motif_logo(
        self, method: str = "probability", export: bool = False, **kwargs
    ):
        """
        Plot the cluster logo.

        Parameters
        ----------
        method: str,  "probability" or "information", default = "probability"
            Create the logo based on either the probability of each AA in the cluster,
            or the information content (shannon surprise) compared against a synthetic
            repertoire.
        export: bool, default = False
            Returns the plot as a BytesIO file-like png object.
        **kwargs:
            Extra arguments passed to the logomaker `Logo`.
        """
        if method not in ["probability", "information"]:
            raise ValueError(
                f'{method} is not a valid method, choose "probability" or "information"'
            )

        input_matrix = self.profile_matrix().transpose().reset_index(drop=True)

        if method == "information":
            bg_matrix = BG_POSITION_MATRICES[self.sequence_len]
            input_matrix = transform_matrix(
                input_matrix,
                from_type="probability",
                to_type="information",
                background=bg_matrix,
            )

        height = kwargs.pop("height", 1)
        width_per_col = kwargs.pop("width_per_col", 0.5)
        color_scheme = kwargs.pop("color_scheme", "NajafabadiEtAl2017")
        highlight_pos = kwargs.pop("highligh_position", None)

        if "ax" not in kwargs:
            fig, ax = plt.subplots(figsize=(width_per_col * self.sequence_len, height))
            ax_in_kwargs = False
        else:
            ax = kwargs.pop("ax")
            ax_in_kwargs = True

        Logo(
            df=input_matrix,
            color_scheme=color_scheme,
            font_name="Arial",
            ax=ax,
            **kwargs,
        )

        plt.axis("off")
        if export:
            b = io.BytesIO()
            plt.savefig(b, format="png")
            plt.close()
            return b
        else:
            plt.close()
            if not ax_in_kwargs:
                return fig
            else:
                return ax

    def kmerize(self, k) -> List[str]:
        return [
            a for b in [list(kmer_iterator(s, k)) for s in self.sequences] for a in b
        ]


class ClusterRepertoire:
    """
    Class for storing a set of clusters (cluster repertoire), calculate their
    useful properties and query them for similarities.
    """

    def __init__(
        self, clusters: List[Cluster], k: Union[int, range, list, str] = None
    ) -> None:
        if isinstance(clusters, pd.Series):
            clusters = clusters.to_numpy()
        self.clusters = clusters

        if k:
            self._k = k
            for c in self.clusters:
                c.k = k
        else:
            _assert_k_equality(clusters)
            self._k = clusters[0].k
        self._k_list = [a for b in [c._k_list for c in self.clusters] for a in b]
        self._k_set = set(self._k_list)

        self._num_perm = 128
        self._num_part = 32
        self._minhash_weights = (0.2, 0.8)
        # immutable:
        self.threshold = 0.2
        self.l = 8

    def __repr__(self) -> str:
        return f"Repertoire ({len(self.clusters)} clusters)"

    def __len__(self) -> int:
        return len(self.clusters)

    def __getitem__(self, i):
        return self.clusters[i]

    def __add__(self, other):  # add two ClusterRepertoires
        new_clst = self.clusters + other.clusters
        idlist = [c.xid for c in new_clst]
        assert len(idlist) == len(set(idlist))
        return ClusterRepertoire(clusters=new_clst)

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def __hash__(self) -> int:
        return hash(self.__repr__())

    def get(self, index: str, default=None) -> Cluster:
        """
        Get cluster from id.
        """
        return self.cluster_dict.get(index, default)

    def __iter__(self):
        return (c for c in self.clusters)

    @classmethod
    def from_clustcr_result(cls, cluster_result: ClusteringResult, **kwargs):
        """
        Alternative constructor to create a ClusterRepertoire from the ClusteringResult
        object generated by ClusTCR.
        """
        cluster_list = [Cluster(seqs) for seqs in cluster_result.cluster_contents()]
        k = kwargs.pop("k", 5)
        return cls(clusters=cluster_list, k=k, **kwargs)

    @property
    def k(self):
        """
        Value used to split cluster sequences into k-mers. Setting k here will set ks for all clusters.

        Additionally; a list or range of values, or 'all' can be passed to use multiple k values.
        """
        _assert_k_equality(self.clusters, self._k)
        return self._k

    @k.setter
    def k(self, k: Union[int, range, list, str]):
        if k != self._k:
            for c in self.clusters:
                c.k = k
            self._k = k
            self._k_list = [a for b in [c._k_list for c in self.clusters] for a in b]
            self._k_set = set(self._k_list)
            self._reset_cached_properties()

    @property
    def num_perm(self):
        return self._num_perm

    @num_perm.setter
    def num_perm(self, num_perm):
        self._num_perm = num_perm
        for c in self.clusters:
            c.num_perm = num_perm
        self._reset_cached_properties()

    @property
    def num_part(self):
        return self._num_part

    @num_part.setter
    def num_part(self, num_part):
        self._num_part = num_part
        self._reset_cached_properties()

    @property
    def minhash_weights(self):
        return self._minhash_weights

    @minhash_weights.setter
    def minhash_weights(self, wts):
        self._minhash_weights = wts
        self._reset_cached_properties()

    def _reset_cached_properties(self):
        for cached_prop in [
            "_kmer_containing_cluster_dict",
            "lsh_index",
            "ensemble_lsh_index",
            "forest_lsh_index",
        ]:
            if cached_prop in self.__dict__:
                del self.__dict__[cached_prop]

    @cached_property
    def cluster_dict(self):
        """
        Dict of shape {cluster_id : Cluster}.
        """
        return {c.xid: c for c in self.clusters}

    @property
    def repertoire_size(self):
        """
        Number of clusters in the ClusterRepertoire.
        """
        return len(self.clusters)

    def as_dataframe(self):
        """
        Return ClusterRepertoire as a dataframe, using the Cluster id as index
        and the Cluster object as a column.
        """
        return pd.DataFrame(
            {"id": [c.xid for c in self.clusters], "cluster_object": self.clusters}
        ).set_index("id", drop=True)

    @property
    def average_pgens(self) -> pd.Series:
        """
        Computes pd.Series with average generation probability for each cluster.

        Note: multiprocessing is used to speed up calculations for large repertoires,
        individual pgen values for the CDR3 sequences are cached and can be retrieved
        as a cluster property.
        """
        if self.repertoire_size < 100:
            res = [c.average_pgen for c in self.clusters]
        else:
            with Pool(4) as p:
                res = p.map(self._pgen_helper, self.clusters)
        return pd.Series(res, index=[c.xid for c in self.clusters])

    def _pgen_helper(self, c: Cluster):
        return c.average_pgen

    @cached_property
    def _kmer_containing_cluster_dict(self):
        """
        Dict with as key all unique kmers in the repertoire, and as values the number of
        clusters in which they appear. Used for tf-idf weighting.
        """
        cnt = Counter([a for b in [list(set(c._k_list)) for c in self] for a in b])
        return dict(cnt)

    @cached_property
    def lsh_index(self) -> MinHashLSH:
        """
        Indexes the ClusterRepertoires using locality sensitive hashing for quick querying
        of similar clusters using jaccard similarity. The resulting index is cached for later use.
        """
        lsh = MinHashLSH(
            threshold=self.threshold,
            num_perm=self.num_perm,
            weights=self.minhash_weights,
        )
        for c in self.clusters:
            lsh.insert(str(c.xid), c.mhash)
        return lsh

    @cached_property
    def ensemble_lsh_index(self) -> MinHashLSHEnsemble:
        """
        A locality sensitive hashing implementation for containment searches (as opposed to jaccard
        similarity). The resulting index is cached for later use.
        """
        lshe = MinHashLSHEnsemble(
            threshold=self.threshold,
            num_perm=self.num_perm,
            num_part=self.num_part,
            weights=self.minhash_weights,
        )
        lshe.index([(str(c.xid), c.mhash, len(set(c._k_list))) for c in self.clusters])
        return lshe

    @cached_property
    def forest_lsh_index(self):
        """
        Indexes using a locality sensitive hashing implementation optimized to query the
        ClusterRepertoire for the top-n clusters based on jaccard similarity. The result
        is cached for later use.
        """
        lshf = MinHashLSHForest(
            num_perm=self.num_perm, l=self.l, weights=self.minhash_weights
        )
        for c in self.clusters:
            lshf.add(c.xid, c.mhash)
        lshf.index()
        return lshf

    def query_similar_clusters(
        self,
        cluster: Cluster,
        method: str = "jaccard",
        n: int = 20,
        return_id: bool = False,
    ) -> List[Cluster]:
        """
        Uses an efficient minhash implementation to quickly query the cluster repertoire
        for clusters similar to the provided cluster.

        Parameters
        ----------
        cluster : Cluster
            The query cluster.
        method : str["jaccard", "containment", "top-n"], default="jaccard"
            Minhash implementation to use. The default "jaccard" optimizes matches
            for maximal jaccard similarity. The containment method normalizes for large
            variations in cluster size. Top-n uses multiple jaccard indexes to return a
            list of size n of top matches.
        n : int, optional, default=20
            The size of return list when using the top-n method.
        return_id:  bool, optional, default=False
            If True: return the id of matching clusters, as opposed to the cluster objects.

        Returns
        -------
        List[Cluster]
            The list of matching clusters using the given method.
        """
        cluster.k = self.k
        if method == "jaccard":
            id_list = self.lsh_index.query(cluster.mhash)

        elif method == "containment":
            l = len(set(cluster._k_list))
            id_list = list(self.ensemble_lsh_index.query(cluster.mhash, l))

        elif method in ["top", "top-n", "top_n"]:
            id_list = self.forest_lsh_index.query(cluster.mhash, k=n)

        else:
            raise ValueError(f"method '{method}' is not a valid method")

        if return_id:
            return [int(xid) for xid in id_list]
        else:
            return [self.cluster_dict.get(int(xid)) for xid in id_list]

    def set_read_counts(self, reads_counter: Union[Counter, dict]) -> None:
        """
        set all cluster read count properties using a dict of shape
        {CDR3: count}
        """
        for c in self.clusters:
            c.read_count = reads_counter


def _assert_k_equality(clusterlist=None, k=None, c: Cluster = None) -> None:
    if clusterlist:
        ks = set([c.k for c in clusterlist])
        if len(ks) != 1:
            raise ValueError(f"All clusters must have same k-value, {ks} were found.")
        if k:
            if k not in ks:
                raise ValueError(f"Cluster k and ClusterRepertoire k do not match")
    if c and k:
        if c.k != k:
            raise ValueError(f"Query k {c.k} does not match ClusterRepertoire k {k}")
