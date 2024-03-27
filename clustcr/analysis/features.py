import numpy as np
import pandas as pd

from clustcr.chem_properties import PHYSCHEM
from .tools import profile_matrix, motif_from_profile, format_chain, get_olga_model

class FeatureGenerator:
    """
    Compute feature matrix that provides a numerical representation of the
    sequences within a cluster. The feature matrix can be used for exploration
    of clustering results, and for downstream machine learning applications.
    
    Calculated features include:
        - entropy
        - size
        - length of CDR3
        - physicochemical features
        - generation probability
    """
    def __init__(self, nodes):
        self.nodes = nodes
        self.clusterids = nodes['cluster'].unique()



    def _calc_variation(self, correction="log"):
        """
        Correction factors:
            - log: 1/log2(n) correction (better for smaller clusters)
            - ssc: small-sample correction (typically used in motif logo construction)
        """
        
        # Correction factors
        cfactors = ["log", "ssc"]
        assert correction in cfactors, "Unknown correction factor '{}', please choose one of the following: {}.".format(correction, cfactors)
        
        # Results
        res = {"cluster":[], "h":[], "size":[], "length":[]}
        
        # Calculate average information content per amino acid position in cluster
        for clust in self.nodes["cluster"].unique():
            
            sequences = self.nodes[self.nodes["cluster"]==clust]["junction_aa"].tolist() # sequences of one cluster
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
            
            res["cluster"].append(clust)
            res["h"].append(np.average(ic))
            res["size"].append(n)
            res["length"].append(l)
        
        return pd.DataFrame(res).set_index("cluster", drop=True)


    
    def _calc_physchem(self):
        """
        Calculate the average physicochemical properties for a CDR3 amino acid sequence.
        Takes a nodelist as input. This can be calculated with the network_clustering() function.
        
        To add a physicochemical property, add a dictionary containing the values for each AA,
        also add physicochemical property to the "physchem_properties" dictionary.
        """
        physchem_properties = PHYSCHEM

        properties = {}
        for seq in self.nodes["junction_aa"]:
            for prop in physchem_properties:
                if prop not in properties:
                    properties[prop] = []
                properties[prop].append(np.average([physchem_properties[prop][aa] for aa in seq]))
        for prop in properties:
            self.nodes[prop] = properties[prop]
        
        cols_1 = [prop + "_avg" for prop in list(physchem_properties.keys())]
        cols_2 = [prop + "_var" for prop in list(physchem_properties.keys())]
        physchemprop = pd.concat([pd.DataFrame([self.nodes.groupby("cluster")[prop].mean() for prop in physchem_properties], index=cols_1).T,
                                  pd.DataFrame([self.nodes.groupby("cluster")[prop].var() for prop in physchem_properties], index=cols_2).T],
                                 axis= 1)
        
        return physchemprop 
    
    
    
    def _calc_pgen(self, chain='B'):
        """
        Calculate the average generation probability of a cluster.
        Pgen calculations are based on the OLGA module.
        """
        # Select alpha or beta Pgen model
        pgen_model = get_olga_model(chain=chain)
        # Compute Pgen
        p = [pgen_model.compute_aa_CDR3_pgen(seq) for seq in self.nodes["junction_aa"]]
        # Format results
        self.nodes["pgen"] = p
        pgenvals = pd.concat([self.nodes.groupby("cluster")["pgen"].mean().rename("pgen_avg"),
                              self.nodes.groupby("cluster")["pgen"].var().rename("pgen_var")],
                             axis=1)
        return pgenvals
        
        
    
    def _combine(self, *args):
        """
        Combine cluster features into one pd.DataFrame.
        """
        return pd.concat([*args], axis=1)
    


    def get_features(self, chain='B', compute_pgen=True):
        """
        Compute feature matrix. Calculating pgen values can be time consuming.
        Therefore, the user has the option to skip this step, by turning the
        compute_pgen parameter to False.
        """
        chain = format_chain(chain)
        aavar = self._calc_variation()
        pchem = self._calc_physchem()
        if compute_pgen:
            pgen = self._calc_pgen(chain=chain)
            return self._combine(aavar, pchem, pgen)
        else:
            return self._combine(aavar, pchem)
    
    
    
    def clustermotif(self, method='standard', cutoff=.7):
        """
        Calculate a consensus motif representation for a set of sequence
        based on the profile matrix.
        """
        clustermotifs = dict()
        for i in self.clusterids:
            sequences = self.nodes[self.nodes['cluster'] == i]['junction_aa'].tolist()
            profile = profile_matrix(sequences)
            motif = motif_from_profile(profile, method=method, cutoff=cutoff)
            clustermotifs[i] = motif
        
        return clustermotifs
