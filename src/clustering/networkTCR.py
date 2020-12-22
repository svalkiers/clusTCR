# -*- coding: utf-8 -*-
"""
author: Sebastiaan Valkiers
"""

import numpy as np
import pandas as pd
import networkx as nx
import markov_clustering as mcl

import modules.olga.load_model as load_model
import modules.olga.generation_probability as pgen

from modules.faiss_clustering import FaissClustering, DistancePairs
from clustering.tools import create_edgelist, profile_matrix, motif_from_profile
from clustering.amino_acid import PHYSCHEM
from load_files.datasets import vdj_small_cdr3


def test_data():
    return vdj_small_cdr3()
    

class Clustering:
    
    def __init__(self, cdr3, method):
        self.cdr3 = cdr3
        self.method = method.upper()
        
        available = ["MCL", 
                     "FAISS",
                     "TWO-STEP"]
        
        assert self.method in available, "Method not available, please choose one of the following methods: \n {}".format(available)
    


    def MCL(self, edgelist = None, mcl_hyper=[1.2,2], outfile=None):
    
        '''
        Perform clustering on a network of CDR3 amino acid sequences with a known hamming distance,
        using the Markov clustering (MCL) algorithm. For more info about the inflation and expansion
        parameters, visit: https://micans.org/mcl/
        '''
        
        if edgelist is None:
            edgelist = create_edgelist(self.cdr3)
        else:
            pass
        
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
    
    
    
    def FAISS(self, cluster_size = 10):
        
        nodelist = {"CDR3":[], "cluster":[]}
        clusters = FaissClustering.cluster(self.cdr3, avg_items_per_cluster = cluster_size)
        for cluster in clusters:
            nodelist["CDR3"].append(cluster)
            nodelist["CDR3"].append(len(nodelist["cluster"]))
            
        return nodelist
    
    
    def TWOSTEP(self, size_of_preclusters = 500):
        '''
        Two-step clustering procedure that combines the speed of the faiss method
        with the accuracy of MCL.
        '''
        
        # Pre-sorting sequences using faiss
        preclust = FaissClustering.cluster(self.cdr3, avg_items_per_cluster = size_of_preclusters)
        
        # Actual clustering using MCL
        initiate = True
        for c in preclust.get_cluster_contents():
            try:
                edges = create_edgelist(c)
                if initiate:
                    nodelist = self.MCL(edges)
                    initiate = False
                else:
                    nodes = self.MCL(edges)
                    nodes["cluster"] = nodes["cluster"] + nodelist["cluster"].max() + 1
                    nodelist = nodelist.append(nodes)
            # If no edges can be found, leave cluster as is
            except nx.NetworkXError:
                cluster = pd.DataFrame({"CDR3" : c,
                                        "cluster" : [nodelist["cluster"].max() + 1] * len(c)})
                nodelist = nodelist.append(cluster)
        
        return nodelist
    


    def TWOSTEP_V2(self,
                   size_of_preclusters = 500,
                   edge_method = 'distance matrix',
                   max_allowed_distance = 5,
                   weighting_scheme = 2):
        '''
        UNDER CONSTRUCTION
        '''
        # Pre-sorting sequences using faiss
        preclust = FaissClustering.cluster(self.cdr3, avg_items_per_cluster = size_of_preclusters)
        
        if edge_method.upper() == 'HASHING':
            
            # Actual clustering using MCL
            initiate = True
            for c in preclust.get_cluster_contents():
                try:
                    edges = create_edgelist(c)
                    if initiate:
                        nodelist = self.MCL(edges)
                        initiate = False
                    else:
                        nodes = self.MCL(edges)
                        nodes["cluster"] = nodes["cluster"] + nodelist["cluster"].max() + 1
                        nodelist = nodelist.append(nodes)
                # If no edges can be found, leave cluster as is
                except nx.NetworkXError:
                    cluster = pd.DataFrame({"CDR3" : c,
                                            "cluster" : [nodelist["cluster"].max() + 1] * len(c)})
                    nodelist = nodelist.append(cluster)
        
        elif edge_method.upper() == 'DISTANCE MATRIX':
            
            distances = DistancePairs.generate(preclust, self.cdr3)
            weighted_edges = distances.get_weighted_edges()
        
        return weighted_edges
    
    
    
    def TCR_clustering(self):
        
        if self.method == 'MCL':
            nodelist = self.MCL()
            
        elif self.method == 'FAISS':
            nodelist = self.FAISS()
            
        elif self.method == 'TWO-STEP':
            nodelist = self.TWOSTEP()
            
        return nodelist
        
            

class Features:
    # Calculate features for clusters
    def __init__(self, nodes):
        self.nodes = nodes
        self.clusterids = nodes['cluster'].unique()



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
        physchem_properties = PHYSCHEM

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
            
            params_file_name = 'modules/olga/default_models/human_T_beta/model_params.txt'
            marginals_file_name = 'modules/olga/default_models/human_T_beta/model_marginals.txt'
            V_anchor_pos_file ='modules/olga/default_models/human_T_beta/V_gene_CDR3_anchors.csv'
            J_anchor_pos_file = 'modules/olga/default_models/human_T_beta/J_gene_CDR3_anchors.csv'
            
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
    
    
    
    def clustermotif(self):
        
        clustermotifs = dict()
        for i in self.clusterids:
            sequences = self.nodes[self.nodes['cluster'] == i]['CDR3'].tolist()
            profile = profile_matrix(sequences)
            motif = motif_from_profile(profile)
            clustermotifs[i] = motif
        
        return clustermotifs
    

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
        conf_mat_t = pd.pivot_table(self.gt,values='count',
                                    index=self.gt["Epitope"],
                                    columns=self.gt["cluster"],
                                    aggfunc=np.sum,
                                    fill_value=0)
        conf_mat_b = pd.pivot_table(self.gt_baseline,
                                    values='count',
                                    index=self.gt_baseline["Epitope"],
                                    columns=self.gt_baseline["cluster"],
                                    aggfunc=np.sum,
                                    fill_value=0)
        
        return conf_mat_t, conf_mat_b
    
    
             
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