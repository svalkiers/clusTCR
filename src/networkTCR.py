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


class clustering:
    
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
        
        assert type(self._set) == set, "Method createNetwork() requires a set."
        
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
        - Write output to file (optional).
        
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
        nodelist = {"CDR3":[], "cluster":[]}
        for c in cluster_ids:
            for m in cluster_ids[c]:
                nodelist["CDR3"].append(m)
                nodelist["cluster"].append(c)
        nodelist = pd.DataFrame(data=nodelist)
        
        # Write to file
        if outfile is not None:
            nodelist.to_csv(outfile, sep="\t", index=False)
            
        return nodelist
    


    def calc_profile(self, sequences):
        '''
        Calculates the profile matrix for a set of sequences (i.e. all cluster members).
        NOTE: this version does not take into account the expected frequency of each amino acid at each position.
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
        profile = pd.DataFrame(profile,index=colnames).T # indices will be columns, because the df is transposed
        
        return profile
    
    
    
    def clustermotif(self, profile):
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
    
    

class metrics:
    # NOTE: You can only calculate these cluster metrics if the CDR3 sequences are labelled (i.e. epitope specificity is known)!
    def __init__(self, nodelist, epidata):
        self.nodelist = nodelist # pd.DataFrame with columns ["CDR3", "cluster"]
        self.epidata = epidata # 'ground truth', pd.DataFrame of CDR3 sequences with corresponding epitope specificity (columns=["CDR3", "Epitope"])
        
        
        
    def retention(self):
        '''
        Cluster retention is the fraction of sequences that has been assigned to any cluster.
        '''
        return len(self.nodelist) / len(self.epidata["CDR3"].unique())
    
        
    
    def purity(self, weighted=True):
        '''
        Cluster purity is fraction of CDR3 sequences within a single cluster that target the same epitope.
        This metric describes the purity of that cluster, in terms of epitope specificity.
        
        This function also provides a baseline estimation of cluster purity, by permuting the cluster assignment
        column and calculating purity of the permuted data using the same procedure.
        '''
        
        # Ensure all values correspond to CDR3s in nodelist and no duplicates remain
        self.epidata = self.epidata[self.epidata["CDR3"].isin(self.nodelist["CDR3"])]
        self.epidata.drop_duplicates(inplace=True)
        
        # Make new column "permuted", which forms a baseline (as if the clustering was random)
        self.nodelist["permuted"] = np.random.permutation(self.nodelist["cluster"])
                
        # Construct joint pd.DataFrame that stores information about cluster and epitope association of CDR3s
        gt = pd.merge(left=self.epidata, right=self.nodelist, on="CDR3")
        
        # Calculate purity of each individual cluster and take average
        # Regular cluster assignments
        pty_reg = [gt[gt["cluster"]==i]["Epitope"].value_counts()[0] / len(gt[gt["cluster"]==i]["CDR3"].unique()) for i in sorted(gt["cluster"].unique())]
        pty_avg_reg = np.average(pty_reg)
        
        # Permuted cluster assignments
        pty_base = [gt[gt["permuted"]==i]["Epitope"].value_counts()[0] / len(gt[gt["permuted"]==i]["CDR3"].unique()) for i in sorted(gt["permuted"].unique())]
        pty_avg_base = np.average(pty_base)
        
        # Recalculate purity by weighting clusters based on their size
        if weighted == True:
            
            # Regular
            size_reg = [len(gt[gt["cluster"]==i]["CDR3"].unique()) for i in sorted(gt["cluster"].unique())]
            pty_reg = [pty_reg[n] * np.log(size_reg[n]) for n in range(len(size_reg))]
            pty_avg_reg = sum(pty_reg) / sum(np.log(size_reg))
            
            # Permuted
            size_per = [len(gt[gt["permuted"]==i]["CDR3"].unique()) for i in sorted(gt["permuted"].unique())]
            pty_base = [pty_base[n] * np.log(size_per[n]) for n in range(len(size_per))]
            pty_avg_base = sum(pty_base) / sum(np.log(size_per))            
        
        return {"Regular":pty_avg_reg, "Baseline":pty_avg_base}