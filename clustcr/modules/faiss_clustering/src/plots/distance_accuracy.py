import matplotlib.pyplot as plt
import numpy as np
from ..distance import make_pairs_with_distance, make_pairs_with_distance_after_sorting


def plot_accuracy(data, dataset_name):
    plt.figure()
    total_pairs_of_distance1 = len(createNetwork(data))

    x = []
    y_sorted_faiss = []
    y_faiss = []

    for items_per_cluster in np.arange(10, 110, 30):
        x.append(items_per_cluster)

        pairs_sorted_faiss, contents = make_pairs_with_distance_after_sorting(data, int(items_per_cluster))
        accuracy = distance1_accuracy(pairs_sorted_faiss, data, total_pairs_of_distance1)
        y_sorted_faiss.append(accuracy)

        pairs_faiss, contents = make_pairs_with_distance(data, int(items_per_cluster))
        accuracy = distance1_accuracy(pairs_faiss, data, total_pairs_of_distance1)
        y_faiss.append(accuracy)

        # Progress
        print(items_per_cluster)

    plt.plot(x, y_faiss, label='Faiss method')
    plt.plot(x, y_sorted_faiss, label='Faiss method (after sorting on length)')

    plt.xlabel('Average amount of items per cluster')
    plt.ylabel('Percentage of pairs found')
    plt.ylim(0, 1)
    plt.legend(loc='best')
    plt.title(f'Finding pairs with distance 1 (in {dataset_name})')
    plt.show()


def distance1_accuracy(pairs, data, total_pairs_of_distance1=None):
    pairs_of_distance1 = [p for p in pairs if p[0] == 1]
    if total_pairs_of_distance1 is None:
        total_pairs_of_distance1 = len(createNetwork(data))
    return len(pairs_of_distance1) / total_pairs_of_distance1


def createNetwork(_set, dist=1, filename=None):
    '''
    Creates a network where nodes are represented by CDR3 sequences and edges are the edit distance (dist) between them.
    The algorithm finds matches by hashing the sequences. This provides accurate results for dist = 1, but is not fully
    accurate for dist > 1.
    '''
    # Hashing
    cdr3hash = dict()
    for cdr in _set:
        for hash in (cdr[::2], cdr[1::2]):
            if hash not in cdr3hash:
                cdr3hash[hash] = set()
            cdr3hash[hash].add(cdr)
    # Generate network
    edgelist = set()
    for hash in cdr3hash:
        if len(cdr3hash[hash]) >= 1:
            for cdr1 in cdr3hash[hash]:
                for cdr2 in cdr3hash[hash]:
                    if cdr1 != cdr2:
                        if cdr1 <= cdr2:
                            if sum(ch1 != ch2 for ch1, ch2 in zip(cdr1, cdr2)) <= dist:
                                edgelist.add(cdr1 + "\t" + cdr2)
    # Write results to file
    if filename is not None:
        with open(filename, 'w') as f:
            for item in edgelist:
                f.write("%s\n" % item)
    return edgelist
