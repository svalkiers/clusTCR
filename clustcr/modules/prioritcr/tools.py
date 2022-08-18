from functools import partial
import time
from typing import Callable, Iterable, Tuple, Union, List, Iterator
from collections import defaultdict
from itertools import chain, repeat, combinations
from multiprocessing import Pool

import numpy as np
from scipy.spatial.distance import squareform
import mmh3


def timed(myfunc):
    # Decorator to keep track of time required to run a function
    def timed(*args, **kwargs):
        start = time.time()
        result = myfunc(*args, **kwargs)
        end = time.time()
        print(f"Total time to run '{myfunc.__name__}': {(end-start):.3f}s")
        return result

    return timed


def kmer_iterator(s: str, k: Union[int, list]) -> Iterator[str]:
    """
    yields all k-mers of lenght k in a sequence s

    additionally, a list or range of ks can be passed

    Examples
    --------
    >>> list(kmer_iterator("ABCDE", 2))
    ['AB', 'BC', 'CD', 'DE']

    >>> list(kmer_iterator("ABCDE", range(2,4)))
    ['AB', 'BC', 'CD', 'DE', 'ABC', 'BCD', 'CDE']

    >>> list(kmer_iterator("ABC", "all"))
    ['A', 'B', 'C', 'AB', 'BC', 'ABC']
    """
    if isinstance(k, int):
        for i in range(len(s) - k + 1):
            yield s[i : i + k]

    elif isinstance(k, (range, list)):
        for k_ in k:
            for i in range(len(s) - k_ + 1):
                yield s[i : i + k_]

    elif k in ["all", "max"]:
        kmax = len(s) + 1
        for k_ in range(1, kmax):
            for i in range(len(s) - k_ + 1):
                yield s[i : i + k_]


def hamming_intersection(set1: set, set2: set):
    """
    returns the intersection of sequences from set1 and set2,
    including sequences in set2 that are a hamming distance of 1 away
    from a sequence in set1
    """
    d = defaultdict(lambda: defaultdict(set))
    # dict(sequence_hash : dict(1 : set[sequences], 2 : set[sequences]))
    for s, i in chain(zip(set1, repeat(1)), zip(set2, repeat(2))):
        for hash in (s[::2], s[1::2]):
            d[hash][i].add(s)

    return {
        s2
        for m in d.values()
        for s1 in m[1]
        for s2 in m[2]
        if len(s1) == len(s2)
        if sum([aa1 != aa2 for aa1, aa2 in zip(s1, s2)]) <= 1
    }
    # iterates through sequences with the same hash, creates set of sequences from s2 if a sequence with hamming distance <= 1 is in the hash list


def hash_kmer(kmer, m):
    """
    Uses MurMurHash for fast hashing of strings
    """
    if m == 32:
        return mmh3.hash(kmer, signed=False)
    elif m == 64:
        return mmh3.hash64(kmer, signed=False)[0]
    elif m == 128:
        return mmh3.hash128(kmer, signed=False)
    else:
        raise ValueError(f"{m} not valid, choose 32,64,128")


def vec_bin_array(arr, m):
    """
    Create binary np array from corresponding int arr
    """
    convert_int = lambda n: [
        int(i) for i in n.to_bytes(m // 8, byteorder="big", signed=False)
    ]
    byte_arr = np.array([convert_int(n) for n in list(arr)], dtype="uint8")
    bit_arr = np.unpackbits(byte_arr).reshape(byte_arr.shape[0], byte_arr.shape[1] * 8)
    return bit_arr


def create_distance_matrix(cr, metric_function: Callable, n_cpus: int = 1):
    """
    Exhaustively create a distance matrix using a given function and cluster
    repertoire.
    """

    if metric_function.__name__ in ["tfidf_cosine", "tifidf_kmer_jaccard"]:
        mf = partial(metric_function, background_a=cr, background_b=cr)
    else:
        mf = metric_function

    if n_cpus > 1:
        with Pool(n_cpus) as pool:
            distance_list = pool.starmap(mf, list(combinations(cr, 2)))
    else:
        distance_list = [mf(a, b) for a, b in combinations(cr, 2)]

    distance_matrix = squareform(1 - np.array(distance_list))
    return distance_matrix


def encode_permutation(arr: np.ndarray, combination: Iterable, partition_size: int):
    """
    For multi-hashtable indexing
    """
    p = partition_size
    if isinstance(combination, int):
        i = combination
        res = np.apply_along_axis(bin2int, 1, arr[:, i * p : i * p + p])
    else:
        res = np.apply_along_axis(
            bin2int, 1, np.hstack([arr[:, i * p : i * p + p] for i in combination])
        )
    return res


def encode_1d_permutation(arr: np.ndarray, combination: tuple, partition_size: int):
    """
    Encode a single permutation
    """
    p = partition_size
    if isinstance(combination, int):
        i = combination
        res = bin2int(arr[i * p : i * p + p])
    else:
        res = bin2int(np.concatenate([arr[i * p : i * p + p] for i in combination]))
    return res


def bin2int(x: np.ndarray):
    """
    Convert binary array to int
    """
    y = 0
    for i, j in enumerate(x):
        y += j << i
    return y


def construct_simhash_index(
    xids: list, hashes=list, permutations=range(4), partition_size=8
):
    """
    Construct multihash table using given permutation scheme
    """
    hashmap = defaultdict(set)
    hashes = np.array(hashes)

    for perm in permutations:
        for xid, hashed_perm in zip(
            xids, encode_permutation(hashes, perm, partition_size)
        ):
            hashmap[hashed_perm].add(xid)
    return hashmap
