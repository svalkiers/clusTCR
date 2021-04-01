def create_edgelist(cdr3, cutoff=2, filename=None):
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

    # Save edgelist to file
    if filename is not None:
        with open(filename, 'w') as f:
            for edge in edgelist:
                f.write('%s\n' % edge)

    return edgelist
