def make_pairs_with_distance(data, cluster_contents):
    pairs = []
    for cluster in cluster_contents:
        cluster_pairs = make_pairs_with_distance_in_cluster(data, cluster)
        pairs.extend(cluster_pairs)
    return pairs


def make_pairs_with_distance_in_cluster(data, cluster):
    pairs = []
    for index, id1 in enumerate(cluster):
        for index2 in range(index + 1, len(cluster)):
            id2 = cluster[index2]
            distance = calculate_distance(data[id1], data[id2])
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

