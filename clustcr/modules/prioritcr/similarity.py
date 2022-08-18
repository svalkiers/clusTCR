from collections import Counter
from math import log
from typing import Callable, Tuple, List
from itertools import repeat, product
from functools import partial

import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.spatial.distance import cosine
from scipy.stats import pearsonr
from Levenshtein import median, median_improve, ratio
from .analysis import Cluster, ClusterRepertoire


def levenshtein_similarity(a: Cluster, b: Cluster) -> float:
    """
    Calculate levenshtein similarity
    """
    centre_a = median_improve(median(a.sequences), a.sequences)
    centre_b = median_improve(median(b.sequences), b.sequences)
    return ratio(centre_a, centre_b)


def length_based(a: Cluster, b: Cluster) -> float:
    """
    Calculates similarity based on relative CDR3 length difference
    """
    l1, l2 = a.sequence_len, b.sequence_len
    return min([l1, l2]) / max([l1, l2])


def sequence_jaccard(a: Cluster, b: Cluster) -> float:
    """
    Calculates similarity based on the Jaccard similarity coeficcient of the
    CDR3 sequences.
    """
    s1, s2 = set(a), set(b)
    return len(s1 & s2) / len(s1 | s2)


def sequence_overlap_coefficient(a: Cluster, b: Cluster) -> float:
    """
    Calculates the similarity of two clusters based on the Szymkiewicz-Simpson
    coefficient. https://en.wikipedia.org/wiki/Overlap_coefficient
    """
    s1, s2 = set(a), set(b)
    return len(s1 & s2) / min(len(s1), len(s2))


def position_freq_matrix_rho(a: Cluster, b: Cluster) -> float:
    """
    Calculates similarity based on the position frequency matrix of each
    cluster. This method considers the AA frequencies for each position as a
    probability distribution, and calculates their pearson correlation position
    by position. It returns the mean rho over all positions.
    """
    _assert_equal_len(a, b)
    pma = a._create_profile_matrix(return_df=False)
    pmb = b._create_profile_matrix(return_df=False)
    return sum([pearsonr(pma[:, i], pmb[:, i])[0] for i in range(a.sequence_len)])


def kmer_jaccard(a: Cluster, b: Cluster) -> float:
    """
    Calculates the Jaccard similarity coefficient of the two k-mer sets of each
    cluster.
    """
    _assert_equal_k(a, b)
    s1, s2 = set(a._k_list), set(b._k_list)
    return len(s1 & s2) / len(s1 | s2)


def kmer_multiset_jaccard(a: Cluster, b: Cluster) -> float:
    """
    Calculates the Jaccard similarity coefficient of the two k-mer bags
    (multisets) of each cluster. In contrast to `kmer_jaccard`, this method
    takes the frequency of each k-mer into account. Importantly, the maximum
    value is 1/2.
    """
    _assert_equal_k(a, b)
    ms1, ms2 = Counter(a._k_list), Counter(b._k_list)
    return sum((ms1 & ms2).values()) / sum((ms1 + ms2).values())


def kmer_overlap_coefficient(a: Cluster, b: Cluster) -> float:
    """
    Calculates the similarity of two clusters based on the Szymkiewicz-Simpson coefficient.
    https://en.wikipedia.org/wiki/Overlap_coefficient
    """
    _assert_equal_k(a, b)
    s1, s2 = set(a._k_list), set(b._k_list)
    return len(s1 & s2) / min(len(s1), len(s2))


def kmer_multiset_overlap_coefficient(a: Cluster, b: Cluster) -> float:
    """
    Calculates the similarity of two clusters based on the Szymkiewicz-Simpson
    coefficient, taking the cardinality of each kmer into account.
    """
    _assert_equal_k(a, b)
    ms1, ms2 = Counter(a._k_list), Counter(b._k_list)
    bag_intersection = sum(
        [min(ms1.get(k, 0), ms2.get(k, 0)) for k in ms2.keys() & ms1.keys()]
    )
    return bag_intersection / min(len(a._k_list), len(b._k_list))


def kmer_cosine(a: Cluster, b: Cluster) -> float:
    """
    Calculates the cosine similarity of two clusters. Two vectors are defined
    for each cluster with length the of the intersection of the two kmer sets,
    and as values the counts of those respective kmers for each cluster. The
    cosine similarity is returned as 1 - the cosine of the angle between the
    vectors.
    """
    ms1, ms2 = Counter(a._k_list), Counter(b._k_list)
    union = set(ms1) | set(ms2)
    a, b = [ms1.get(kmer, 0) for kmer in union], [ms2.get(kmer, 0) for kmer in union]
    return 1 - cosine(a, b)


def shash32_hamming_sim(a: Cluster, b: Cluster) -> float:
    hd = np.count_nonzero(a.shash_32 != b.shash_32)
    return (32 - hd) / 32


def shash64_hamming_sim(a: Cluster, b: Cluster) -> float:
    hd = np.count_nonzero(a.shash_64 != b.shash_64)
    return (64 - hd) / 64


def shash128_hamming_sim(a: Cluster, b: Cluster) -> float:
    hd = np.count_nonzero(a.shash_128 != b.shash_128)
    return (128 - hd) / 128


def tfidf_cosine(
    a: Cluster,
    b: Cluster,
    background_a: ClusterRepertoire,
    background_b: ClusterRepertoire = None,
) -> float:
    """
    Calculates the similarity of two clusters, using a tf-idf weighting scheme
    which gives a higher weight to commonly occuring k-mers in the clusters,
    yet that are rare in the repertoire.

    Parameters
    ----------
    a: Cluster
    b: Cluster
        Two clusters of which similarity will be calculated.
    background_a : ClusterRepertoire
        The cluster repertoire used to calculate weightings for cluster a.
    background_b : ClusterRepertoire, optional
        The cluster repertoire used to calculate weightings for cluster b. If no repertoire is passed,
        `background_a` will be used instead.
    """
    _assert_equal_k(a, b)

    # calculate tf-idf weightings
    d1 = _freq_weighted_preprocessing(a, background_a)
    d2 = _freq_weighted_preprocessing(b, background_b if background_b else background_a)

    # compute cosine similarity
    kmers = set(list(d1.keys()) + list(d2.keys()))
    weight_vec_1 = [d1.get(k, 0) for k in kmers]
    weight_vec_2 = [d2.get(k, 0) for k in kmers]
    return 1 - cosine(weight_vec_1, weight_vec_2)


def tfidf_kmer_jaccard(
    a: Cluster,
    b: Cluster,
    background_a: ClusterRepertoire,
    background_b: ClusterRepertoire = None,
) -> float:

    _assert_equal_k(a, b)
    d1 = _freq_weighted_preprocessing(a, background_a)
    d2 = _freq_weighted_preprocessing(b, background_b if background_b else background_a)
    intersection = sum([d1.get(k) + d2.get(k) for k in set(d1) & set(d2)])
    union = sum(d1.values()) + sum(d2.values())
    return intersection / union


def _term_frequency(c: Cluster):
    """
    returns the frequency of terms t in document/cluster D
    """
    total_count = len(c._k_list)
    return {k: v / total_count for k, v in Counter(c._k_list).items()}


def _inverse_document_frequency(k_term: str, bg: ClusterRepertoire):
    """
    calculates log scaled inverse fraction of clusters that
    contain the given k-mer
    """
    N = len(bg)
    n_clusters_containing_term = bg._kmer_containing_cluster_dict.get(k_term, 0)
    return log(N / (1 + n_clusters_containing_term))


def _freq_weighted_preprocessing(c: Cluster, cr: ClusterRepertoire) -> dict:
    """
    calculates kmer tf-idf weighting
    """
    return {
        k: v * _inverse_document_frequency(k, cr) for k, v in _term_frequency(c).items()
    }


def _assert_equal_k(a: Cluster, b: Cluster) -> None:
    if a.k != b.k:
        raise ValueError("both clusters must be initiated with the same k-value")


def _assert_equal_len(a: Cluster, b: Cluster) -> None:
    if a.sequence_len != b.sequence_len:
        raise ValueError("both clusters must have the same CDR3 length")


class ClusterRepertoireSimilarity:
    """
    Class for the quick identification of similar clusters between
    two repertoires, or within the same repertoire.
    """

    def __init__(
        self,
        cluster_repertoire_a: ClusterRepertoire,
        cluster_repertoire_b: ClusterRepertoire = None,
    ) -> None:
        """
        Parameters:
        -----------
        cluster_repertoire_a : ClusterRepertoire
        cluster_repertoire_b : ClusterRepertoire
            The two cluster repertoires of which similar clusters are determined.
            If the cluster_repertoire_b parameter is not passed, or if the same repertoire is passed
            twice, the class will compare clusters within a single repertoire.
        """
        if cluster_repertoire_a == cluster_repertoire_b or cluster_repertoire_b is None:
            self.a, self.b = cluster_repertoire_a, cluster_repertoire_a
        else:
            l = [cluster_repertoire_a, cluster_repertoire_b]
            l.sort(key=len)
            self.a, self.b = l

    def find_matching_clusters(
        self,
        first_pass_method: str = "top_n",
        first_pass_n: int = 26,
        second_pass_similarity_func: Callable = kmer_jaccard,
        second_pass_threshold: float = 0,
        second_pass_top_n: int = None,
        for_network: bool = False,
        remove_identical: bool = True,
    ):
        """
        Uses a two pass approach to quickly retrieve matching clusters. The first step
        reduces the complexity by using locality sensitive hashing functions to query each cluster
        for possible matches. The similarity scores between each cluster and its matches is
        then determined in the second step.

        Parameters
        ----------
        first_pass_method : str, default='top_n'
            Method to use for the first pass. The default 'top-n' uses a MinHash LSH Forest index
            to query each clusters to its most similar matches. Alternatively, use 'jaccard' or 'containment'.
        first_pass_n : int, default = 26
            Number of clusters to query in the first pass. Higher n reduces the chance of false negative matches,
            but significantly increases computation time, especially when using a computationally expensive
            second pass metric.
        second_pass_similarity_func : Callable, default = kmer_jaccard
            A similarity function to use to quantify cluster similarity. See those defined in the
            prioritcr.similarity module
        second_pass_threshold : float, optional
            If passed, only matches with a similarity equal to or higher than this value will be returned.
        second_pass_top_n : int, optional
            If passed, only the top n matches for each cluster will be returned.
        for_network : bool, optional
            If passed, the function will return an edge list using the cluster ids for network visualization.
        remove_identical : bool, default = True
            This removes clusters with an identical id from the result df.

        """
        fpm = self._first_pass_matches(
            method=first_pass_method, n=first_pass_n, remove_identical=remove_identical
        )

        spm = self._second_pass_matches(
            first_pass_matches=fpm,
            similarity_function=second_pass_similarity_func,
            threshold=second_pass_threshold,
        )

        cluster, match = zip(*spm.keys())
        df = pd.DataFrame(
            {"cluster": cluster, "match": match, "similarity": spm.values()}
        )

        if second_pass_top_n:
            df_g = df.groupby("cluster").agg(list)
            df_g[["match", "similarity"]] = df_g.apply(
                lambda r: self._get_n_best_matches(
                    r["match"], r["similarity"], second_pass_top_n
                ),
                axis=1,
                result_type="expand",
            )
            df = df_g.explode(["match", "similarity"]).reset_index()

        if for_network:
            return [(int(a), int(b), 1 - c) for a, b, c in df.to_numpy()]

        c_xid, m_xid, score = df.to_numpy().T
        return pd.DataFrame(
            {
                "cluster": [self.a.get(xid) for xid in c_xid],
                "match": [self.b.get(xid) for xid in m_xid],
                "similarity": score,
            }
        )

    def _get_n_best_matches(self, matches, similarity, n):
        if len(similarity) >= n:
            maxindices = np.argpartition(similarity, -n)[-n:]
            return np.array(matches)[maxindices], np.array(similarity)[maxindices]
        else:
            return matches, similarity

    def _first_pass_matches(
        self,
        method: str = "top_n",
        n: int = 26,
        remove_identical: bool = False,
    ) -> List[Tuple[int, int]]:

        query_func = partial(
            self.b.query_similar_clusters,
            method=method,
            n=n,
            return_id=True,
        )
        matches = [
            zip(repeat(c.xid), query_func(cluster=c))
            for c in tqdm(self.a, desc="First pass matches")
        ]
        matches = [a for b in matches for a in b]
        if remove_identical:
            return [m for m in matches if m[0] != m[1]]
        return matches

    def _second_pass_matches(
        self,
        first_pass_matches: List[Tuple[int]],
        similarity_function: Callable,
        threshold: float,
    ) -> dict:

        # TODO implement multiprocessing here

        if similarity_function.__name__ in ["tfidf_cosine", "tfidf_kmer_jaccard"]:
            similarity_function = partial(
                similarity_function, background_a=self.a, background_b=self.b
            )

        res = {}
        for pair in tqdm(first_pass_matches, desc="Second pass matches"):
            idxa, idxb = pair
            if (idxb, idxa) in res:
                res[pair] = res[(idxb, idxa)]
            else:
                ca, cb = self.a.get(idxa), self.b.get(idxb)
                res[pair] = similarity_function(ca, cb)

        res = {k: v for k, v in res.items() if v >= threshold}
        return res
