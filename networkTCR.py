# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 15:19:50 2020

@author: Sebastiaan Valkiers

Contact: sebastiaan.valkiers@uantwerpen.be
"""

import numpy as np
import pandas as pd
import networkx as nx
import markov_clustering as mcl


class networkTCR:
    
    def __init__(self, _set):
        self._set = _set
        
        

    def createNetwork(self, dist=1, filename=None):
    
        '''
        Creates a network where nodes are represented by CDR3 sequences and edges are the edit distance (dist) between them.
        This is a modified version of the original algorithm written by Pieter Meysman.
        
        Parameters
        ----------
        _set : set, list, array
            Collection of CDR3 sequences.
        dist : int; optional
            Hamming distance between two sequences. The default is 1.
        filename : string, optional
            Name of outfile. The default is None.

        Returns
        -------
        edgelist : set
            Set of CDR3 pairs that have an edit distance = dist.
        '''
        
        # Hashing
        cdr3hash = dict()
        for cdr in self._set:
            for hash in (cdr[::2], cdr[1::2]):
                if hash not in cdr3hash:
                    cdr3hash[hash] = set()
                cdr3hash[hash].add(cdr)
                
        # Generate network
        self.edgelist = set()
        for hash in cdr3hash:
            if len(cdr3hash[hash]) >= 1:
                for cdr1 in cdr3hash[hash]:
                    for cdr2 in cdr3hash[hash]:
                        if cdr1 != cdr2:
                            if cdr1 <= cdr2:
                                if sum(ch1 != ch2 for ch1, ch2 in zip(cdr1, cdr2)) <= dist:
                                    self.edgelist.add(cdr1 + "\t" + cdr2)

        # Write results to file
        if filename is not None:
            with open(filename, 'w') as f:
                for item in self.edgelist:
                    f.write("%s\n" % item)

        return self.edgelist



    def network_clustering(self, edgelist, mcl_hyper=[1.2,2], outfile=None):
    
        '''
        Perform clustering on a network of CDR3 amino acid sequences with a known hamming distance,
        using the Markov clustering (MCL) algorithm. For more info about the inflation and expansion
        parameters, visit: https://micans.org/mcl/
        
        Function performs multiple steps:
        - Generate a nx network from a set of edges
        - Run the MCL algorithm, adjust hyperparameters if necessary (mcl_hyper[inflation,expansion]).
        - Map the cluster numbers back to the input sequences.
        - Update network with attributes (cluster ids).
        - Write output to file.
        
        The output file can be visualized using dedicated network visalisation software such as CytoScape.
        '''
        
        # Generate network using nx
        G = nx.parse_adjlist(edgelist, nodetype=str)
        m = nx.to_scipy_sparse_matrix(G)
        
        # Run MCL
        result = mcl.run_mcl(m, inflation=mcl_hyper[0], expansion=mcl_hyper[1])
        clusters = mcl.get_clusters(result)
        identifiers = list(G.nodes())
        
        # Map cluster ids back to seqs
        cluster_ids = dict()
        for i in range(len(clusters)):
            cluster_ids[i] = list(identifiers[i] for i in clusters[i])
            
        # Generate nodelist
        nodelist = {"TCR":[], "cluster":[]}
        for c in cluster_ids:
            for m in cluster_ids[c]:
                nodelist["TCR"].append(m)
                nodelist["cluster"].append(c)
        nodelist = pd.DataFrame(data=nodelist)
        
        # Write to file
        if outfile is not None:
            nodelist.to_csv(outfile, sep="\t", index=False)
            
        return nodelist
    


    def calc_profile(self, sequences):
        '''
        '''
        assert type(sequences) == list

        # Amino acid alphabet
        alphabet = "ARNDCQEGHILKMFPSTWYV"

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
        profile = pd.DataFrame(profile,index=colnames).T
        
        return profile
    
    
    
    def clustermotif(self, profile):
        '''
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