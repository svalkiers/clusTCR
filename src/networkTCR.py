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

import olga.load_model as load_model
import olga.generation_probability as pgen


class Dataloader:
    # infile variable should provide the path to the benchmarking file (e.g. VDJdb)
    def __init__(self, infile):
        self.infile = infile
        
        
        
    def read_bm_file(self, q=None):
        # NOTE: q-score is only used for VDJdb data.
        bm = pd.read_csv(self.infile, sep="\t")
        if q is not None:
            bm = bm[bm["Score"]>=q]
        
        return bm[["CDR3", "Epitope"]].reset_index(drop=True)



class Clustering:
    
    def __init__(self, _set):
        # Providing a set guarantees no duplicates.
        self._set = _set
        assert type(self._set) == set, "Collection of CDR3 sequences must be a set. Convert input using set()."
        
        

    def create_network(self, dist=1, filename=None):
    
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
    
    

class Features:
    # Calculate features for clusters
    def __init__(self, nodes):
        self.nodes = nodes



    def calc_variation(self, correction="log"):
        '''
        Correction factors:
            - log: 1/log2(n) correction (better for smaller clusters)
            - ssc: small-sample correction (typically used in motif logo construction)
        '''
        
        # Correction factors
        cfactors = ["log", "ssc"]
        assert correction in cfactors, "Unknown correction factor '{}', please choose one of the following: {}.".format(correction, cfactors)
        
        # Results
        res = {"ic":[], "size":[], "length":[]}
        
        # Calculate average information content per amino acid position in cluster
        for clust in self.nodes["cluster"].unique():
            
            sequences = self.nodes[self.nodes["cluster"]==clust]["CDR3"].tolist() # sequences of one cluster
            n = len(sequences) # size of cluster
            l = len(sequences[0][1:-1]) # CDR3 length (ignoring pos 1 and -1)
            ic = [] # information content
            en = (1/np.log(2))*((20-1)/(2*n)) # small-sample correction
            
            # Make sure to proceed only if all sequences in the cluster have equal length
            if all(len(seq) == len(sequences[0]) for seq in sequences) is False:
                
                # On the rare occasion that a cluster contains sequences of inequal length.
                # Typically, there is/are only one (or very few) sequence(s) that differ from the avg. CDR3 length in the cluster.
                # Therefore, we use the length of the highest proportion of sequences as the standard, and delete all others.
                s = []
                for i in sequences:
                    s.append(len(i))
                k = pd.Series(s).value_counts().index[0] # Standard cluster length
                for j in sequences:
                    if len(j) != k:
                        del sequences[sequences.index(j)] # Delete all sequences that differ from k in length.
            
            else:
                
                for pos in range(1,l+1): # iterate over position 1 to -1 (leaving out first C and last F)
                    ic_aa_pos = [] # initiate variable to carry information content per position
                    psc = [seq[pos] for seq in sequences] # get first AA of each sequence at position "pos
                    
                    for aa in pd.Series(psc).value_counts():
                        e = (aa/n)*(np.log2(aa/n)) # individual terms of Shannon entropy
                        ic_aa_pos.append(e)
                        
                    information_content = -sum(ic_aa_pos)
                    if correction == "log":
                        ic.append(information_content / np.log2(n))
                    elif correction == "ssc":
                        ic.append(np.log2(20) - (information_content + en))                        
                    
            res["ic"].append(np.average(ic))
            res["size"].append(n)
            res["length"].append(l)
        
        return pd.DataFrame(res)


    
    def calc_physchem(self):
        '''
        Calculate the average physicochemical properties for a CDR3 amino acid sequence.
        Takes a nodelist as input. This can be calculated with the network_clustering() function.
        
        To add a physicochemical property, add a dictionary containing the values for each AA,
        also add physicochemical property to the "physchem_properties" dictionary.
        '''
        
        basicity = {
        'A': 206.4, 'B': 210.7, 'C': 206.2, 'D': 208.6, 'E': 215.6, 'F': 212.1,
        'G': 202.7, 'H': 223.7, 'I': 210.8, 'K': 221.8, 'L': 209.6, 'M': 213.3,
        'N': 212.8, 'P': 214.4, 'Q': 214.2, 'R': 237.0, 'S': 207.6, 'T': 211.7,
        'V': 208.7, 'W': 216.1, 'X': 210.2, 'Y': 213.1, 'Z': 214.9
        }
    
        hydrophobicity = {
            'A': 0.16, 'B': -3.14, 'C': 2.50, 'D': -2.49, 'E': -1.50, 'F': 5.00,
            'G': -3.31, 'H': -4.63, 'I': 4.41, 'K': -5.00, 'L': 4.76, 'M': 3.23,
            'N': -3.79, 'P': -4.92, 'Q': -2.76, 'R': -2.77, 'S': -2.85, 'T': -1.08,
            'V': 3.02, 'W': 4.88, 'X': 4.59, 'Y': 2.00, 'Z': -2.13
        }
    
        helicity = {
            'A': 1.24, 'B': 0.92, 'C': 0.79, 'D': 0.89, 'E': 0.85, 'F': 1.26,
            'G': 1.15, 'H': 0.97, 'I': 1.29, 'K': 0.88, 'L': 1.28, 'M': 1.22,
            'N': 0.94, 'P': 0.57, 'Q': 0.96, 'R': 0.95, 'S': 1.00, 'T': 1.09,
            'V': 1.27, 'W': 1.07, 'X': 1.29, 'Y': 1.11, 'Z': 0.91
        }
    
        mutation_stability = {
            'A': 13, 'C': 52, 'D': 11, 'E': 12, 'F': 32, 'G': 27, 'H': 15, 'I': 10,
            'K': 24, 'L': 34, 'M': 6, 'N': 6, 'P': 20, 'Q': 10, 'R': 17, 'S': 10,
            'T': 11, 'V': 17, 'W': 55, 'Y': 31
        }

        physchem_properties = {'basicity': basicity,
                               'hydrophobicity': hydrophobicity,
                               'helicity': helicity,
                               'mutation stability': mutation_stability}

        properties = []
        for seq in self.nodes["CDR3"]:
            properties.append([np.average([physchem_properties[prop][aa] for aa in seq]) for prop in physchem_properties])
        self.nodes[list(physchem_properties.keys())] = properties
        
        cols_1 = [prop + "_avg" for prop in list(physchem_properties.keys())]
        cols_2 = [prop + "_var" for prop in list(physchem_properties.keys())]
        physchemprop = pd.concat([pd.DataFrame([self.nodes.groupby("cluster")[prop].mean() for prop in physchem_properties], index=cols_1).T,
                                  pd.DataFrame([self.nodes.groupby("cluster")[prop].mean() for prop in physchem_properties], index=cols_2).T],
                                 axis= 1)
        
        return physchemprop 
    
    
    
    def calc_pgen(self):
        '''
        NOTE: this method requires the Python 3 fork of OLGA (https://github.com/dhmay/OLGA)!
        By default, OLGA is installed within this repo.
        '''
        
            
        print("\nCalculating generation probabilities may take a while. Are you sure you want to continue?")
        user_input = input("Confirm: [Y/N] ")
        
        if user_input.lower() in ("y", "yes"):
            
            params_file_name = 'olga/default_models/human_T_beta/model_params.txt'
            marginals_file_name = 'olga/default_models/human_T_beta/model_marginals.txt'
            V_anchor_pos_file ='olga/default_models/human_T_beta/V_gene_CDR3_anchors.csv'
            J_anchor_pos_file = 'olga/default_models/human_T_beta/J_gene_CDR3_anchors.csv'
            
            genomic_data = load_model.GenomicDataVDJ()
            genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
            
            generative_model = load_model.GenerativeModelVDJ()
            generative_model.load_and_process_igor_model(marginals_file_name)
            
            pgen_model = pgen.GenerationProbabilityVDJ(generative_model, genomic_data)
            
            p = [pgen_model.compute_aa_CDR3_pgen(seq) for seq in self.nodes["CDR3"]]
            self.nodes["pgen"] = p
            
            pgenvals = pd.concat([self.nodes.groupby("cluster")["pgen"].mean().rename("pgen_avg"),
                                  self.nodes.groupby("cluster")["pgen"].var().rename("pgen_var")],
                                 axis=1)
            
            return pgenvals
            
        elif user_input.lower() in ("n", "no"):
            return None
        else:
            print(f"Error: Unrecognised input '{user_input}'.")
            return None
        
        
    
    
    def combine(self, *args):
        return pd.concat([*args], axis=1)
    
    

class Metrics:
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