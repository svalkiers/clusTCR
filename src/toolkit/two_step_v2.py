# -*- coding: utf-8 -*-
"""
author: Sebastiaan Valkiers
"""

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