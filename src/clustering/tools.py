import pandas as pd
import numpy as np

from clustering.amino_acid import AALPHABET


def create_edgelist(cdr3):
    '''
    Create tab-separated edgelist of edges with HD = 1, from a set of sequences.    
    '''
    # Set makes sure there are no dupes
    cdr3 = set(cdr3)
    
    # Hashing
    cdr3hash = dict()
    for cdr in cdr3:
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
                            if sum(ch1 != ch2 for ch1, ch2 in zip(cdr1, cdr2)) <= 1:
                                edgelist.add(cdr1 + "\t" + cdr2)

    return edgelist


def profile_matrix(sequences : list):
    '''
    Calculates the profile matrix for a set of sequences (i.e. all cluster members).
    NOTE: this version does not take into account the expected frequency of each amino acid at each position.
    '''

    # Amino acid alphabet
    alphabet = AALPHABET

    # Initiate profile matrix with zeros
    profile = {}
    for aa in alphabet:
        profile[aa] = [0] * len(sequences[0])
    
    # Fill in profile matrix
    for pos in range(len(sequences[0])):
        psc = pd.Series([seq[pos] for seq in sequences]).value_counts()
        for i in psc.index:
            profile[i][pos] = np.round(psc.loc[i] / len(sequences),2)
    
    # Generate output as a pd.DataFrame
    colnames = ["p" + str(p) for p in range(len(sequences[0]))]        
    profile = pd.DataFrame(profile,index=colnames).T # indices will be columns, because the df is transposed
    
    return profile


def motif_from_profile(profile):
    '''
    Generate consensus sequence motif from a profile matrix.
    Square brackets [...] indicate multiple aa possibilities at that position.
    X represents any aa.
    '''
    
    consensus = ''
    for col in profile.columns:
        if profile[col].max() > .5:
            consensus += profile[col].idxmax()
        elif sum(profile[col].nlargest(2)) >= .5:
            if profile[col].nlargest(2)[0] >= 2 * profile[col].nlargest(2)[1]:
                consensus += profile[col].idxmax()
            else:
                char = "[" + ''.join(profile[col].nlargest(2).index) + "]"
                consensus += char
        else:
            consensus += "X"
    
    return consensus