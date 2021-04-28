import time
import pandas as pd
import numpy as np

from Levenshtein import distance as lev
from clustcr import datasets
from itertools import combinations

AALPHABET = 'ARNDCQEGHILKMFPSTWYV'

def timeit(myfunc):
    # Decorator to keep track of time required to run a function
    def timed(*args, **kwargs):
        start = time.time()
        result = myfunc(*args, **kwargs)
        end = time.time()
        print(f'Total time to run \'{myfunc.__name__}\': {(end-start):.3f}s')
        return result
    return timed

def create_mapping(characters):
    mapping = {}
    for i,j in enumerate(characters):
        mapping[j] = i
    return mapping

def map_str_to_float(matrix, mapping):
    for seq in matrix:
        for pos, letter in enumerate(seq):
            seq[pos] = mapping[letter]
    return np.asarray(matrix)

def hamming_distance(X):
    return (X[:, None, :] != X).sum(2)

@timeit
def hashing_method(cdr3, cutoff=1):
    '''
    Create tab-separated edgelist of edges with HD = 1, from a set of sequences.    
    '''
    # Set makes sure there are no dupes
    cdr3 = set(cdr3)

    # Hashing
    cdr3hash = dict()
    for cdr in cdr3:
        hashes = [cdr[i::cutoff + 1] for i in range(cutoff + 1)]
        for hash in hashes:
            if hash not in cdr3hash:
                cdr3hash[hash] = set()
            cdr3hash[hash].add(cdr)

    # Generate network
    edgelist = set()
    for hash in cdr3hash:
        if len(cdr3hash[hash]) <= 1:
            continue
        for cdr1 in cdr3hash[hash]:
            for cdr2 in cdr3hash[hash]:
                if cdr1 != cdr2 \
                        and cdr1 <= cdr2 \
                        and sum(ch1 != ch2 for ch1, ch2 in zip(cdr1, cdr2)) <= cutoff:
                    edgelist.add(cdr1 + "\t" + cdr2)

    return edgelist

@timeit
def hamming_method(cdr3 : pd.Series, cutoff=1):
    # Initiate edge list
    edgelist = set()   
    for i in cdr3.str.len().unique():
        # Start by sorting sequences according to length
        len_sorted = cdr3[cdr3.str.len()==i]
        # Convert each sequence to a list of characters
        m = [list(sequence) for sequence in list(len_sorted)]
        # Convert letter characters to floats
        m = map_str_to_float(m, create_mapping(AALPHABET))
        # Compute hamming distance matrix
        m = hamming_distance(m)
        # Fill lower half and diagonal of the matrix with nan
        m = m.astype(float)
        m[np.arange(m.shape[0])[:,None] > np.arange(m.shape[1])] = np.nan
        np.fill_diagonal(m, np.nan)
        # Extract indices of sequence pairs with HD == d
        for d in range(cutoff+1):
            indices = np.argwhere(m == d)
            # Add edges to edgelist
            for idx in indices:
                node_1 = len_sorted.iloc[idx[0]]
                node_2 = len_sorted.iloc[idx[1]]
                edge = node_1 + "\t" + node_2 + "\t" + str(d) 
                edgelist.add(edge)
    return edgelist

@timeit
def levenshtein_method(cdr3, cutoff=1):
    edgelist = set()
    combos = [comb for comb in combinations(list(cdr3), 2)]
    for combo in combos:
        d = lev(combo[0],combo[1])
        if d <= cutoff:
            edgelist.add(combo[0] + "\t" + combo[1] + "\t" + str(d))
    return edgelist

if __name__=="__main__":
    data = datasets.vdjdb_beta()
    for i in range(1000, 10001, 1000):
        print(i)
        sample = data.sample(i)
        edges_hash = hashing_method(sample, cutoff=3)
        edges_hamm = hamming_method(sample, cutoff=3)
        edges_ld = levenshtein_method(sample, cutoff=3)
