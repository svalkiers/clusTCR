from src.cluster import profile_cluster
from src.io.output import get_cluster_contents, cluster_contents_to_csv
from src.io.datasets import *

import json


def main():
    cdr3, epi = vdj()
    pairs, content = make_pairs_with_distance_after_sorting(cdr3, 500)
    print('Accuracy:', distance1_accuracy(pairs, cdr3))
    cluster_contents_to_csv(content, cdr3, 'clusters.csv')
    with open('distances.json', 'w') as f:
        json.dump(pairs_to_dict(pairs), f)


def make_pairs_with_distance(data, items_per_cluster=10):
    clustering = profile_cluster(data, items_per_cluster=items_per_cluster)
    cluster_contents = get_cluster_contents(clustering)
    pairs = []
    for cluster in cluster_contents:
        cluster_pairs = make_pairs_with_distance_in_cluster(data, cluster)
        pairs.extend(cluster_pairs)
    return pairs, cluster_contents


def make_pairs_with_distance_after_sorting(cdr3, items_per_cluster):
    df = cdr3.to_frame(name='CDR3')
    df['length'] = df['CDR3'].str.len()
    pairs = []
    cluster_contents = []
    for name, group in df.groupby('length'):
        cdr = group['CDR3']
        new_pairs, new_contents = make_pairs_with_distance(cdr, items_per_cluster)
        pairs.extend(new_pairs)
        cluster_contents.extend(new_contents)
    return pairs, cluster_contents


def pairs_to_dict(pairs, data=None):
    """ Supply data if you want the sequence to be in the dict instead of the id """
    d = {}
    for pair in pairs:
        distance = pair[0]
        if distance not in d:
            d[distance] = []
        seq1, seq2 = int(pair[1]), int(pair[2])
        if data is not None:
            seq1, seq2 = data[seq1], data[seq2]
        d[distance].append((seq1, seq2))
    return d


def make_pairs_with_distance_in_cluster(data, cluster):
    pairs = []
    for index, id1 in enumerate(cluster):
        for index2 in range(index + 1, len(cluster)):
            id2 = cluster[index2]
            seq1 = data[id1]
            seq2 = data[id2]
            distance = calculate_distance(seq1, seq2)
            pairs.append((distance, id1, id2))
    return pairs


def calculate_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        return 1000
    if seq1 is seq2:
        return 0

    differences = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            differences += 1
    return differences


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


if __name__ == '__main__':
    main()
