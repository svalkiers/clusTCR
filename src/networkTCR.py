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
        # Providing a set guarantees no duplicates.
        self._set = _set
        assert type(self._set) == set, "Collection of CDR3 sequences must be a set. Convert input using set()."
        
        

    def createNetwork(self, dist=1, filename=None):
    
        '''
        Creates a network where nodes are represented by CDR3 sequences and edges are the edit distance (dist) between them.
        The algorithm finds matches by hashing the sequences. This provides accurate results for dist = 1, but is not fully
        accurate for dist > 1.
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
        
        The output file can be visualized using dedicated network visalisation software such as Cytoscape.
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
        
        # Ensure all values correspond to CDR3s in nodelist and no duplicates remain
        self.gt = self.epidata[self.epidata["CDR3"].isin(self.nodelist["CDR3"])]
        self.gt.drop_duplicates(inplace=True)
        
        # Construct joint pd.DataFrame that stores information about cluster and epitope association of CDR3s
        self.gt = pd.merge(left=self.epidata, right=self.nodelist, on="CDR3")
        
        # Make a copy and permute cluster assignment column, this provides a baseline for comparison    
        self.gt_baseline = self.gt.copy()
        self.gt_baseline["cluster"] = np.random.permutation(self.gt_baseline["cluster"])
    
    
    
    def calc_confmat(self):
        '''
        Construct confusion matrices for true and baseline.
        '''
        self.gt["count"] = 1
        self.gt_baseline["count"] = 1
        conf_mat_r = pd.pivot_table(self.gt,values='count', index=self.gt["Epitope"], columns=self.gt["cluster"], aggfunc=np.sum, fill_value=0)
        conf_mat_b = pd.pivot_table(self.gt_baseline, values='count', index=self.gt_baseline["Epitope"], columns=self.gt_baseline["cluster"], aggfunc=np.sum, fill_value=0)
        
        return conf_mat_r, conf_mat_b
    
    
             
    def retention(self):
        '''
        Cluster retention is the fraction of sequences that has been assigned to any cluster.
        '''
        return len(self.nodelist) / len(self.epidata["CDR3"].unique())
    
    
    
    def purity(self, conf_mat=None):
        '''
        Method that estimates the precision of the solution.
        We assigned each cluster to the most common epitope.
        All other epitopes in the same cluster are considered false positives.
        '''
        
        if conf_mat is None:
            conf_mat = self.calc_confmat()
        
        hits_t = np.sum(conf_mat[0].apply(np.max,axis=0))
        hits_b = np.sum(conf_mat[1].apply(np.max,axis=0))
        
        return {"True":hits_t/np.sum(conf_mat[0].values,axis=None), "Baseline":hits_b/np.sum(conf_mat[1].values,axis=None)}
    
    
    
    def consistency(self, conf_mat=None):
        '''
        Method that pretends that we solved a supervised problem where each cluster corresponds to a single epitope.
        Returns the accuracy of the best solution.
        '''
        
        if conf_mat is None:
            conf_mat = self.calc_confmat()
        
        #Define recursive function that finds the best fit for the diagonal
        def rec_max(mat):
            high  = mat.max().max()
            col = mat.max().idxmax()
            row = mat[col].idxmax()
            
            if(len(mat.index) > 1 and len(mat.columns) > 1):
                high = high + rec_max(mat.drop(row,0).drop(col,1))
            
            return high
        
        return {"True":rec_max(conf_mat[0])/len(self.gt), "Baseline":rec_max(conf_mat[1])/len(self.gt_baseline)}