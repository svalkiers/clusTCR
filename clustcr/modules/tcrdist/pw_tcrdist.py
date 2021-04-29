import pandas as pd
import numpy as np
import networkx as nx
import os

from tcrdist.repertoire import TCRrep
from sklearn.cluster import DBSCAN, KMeans

from clustcr.clustering.clustering import ClusteringResult
from clustcr.input.vdjdb import parse_vdjdb


def TCRDist(data=None, chain='beta', sparse=False):
    '''
    Calculate distances between TCR sequences, using the TCRDist metric
    described in Dash et al. 2017 Nature.

    Parameters
    ----------
    data : pandas.DataFrame, optional
        Input dataframe containing information about
        CDR3 sequence, V gene and antigen specificity. The default is None.
    chain : str, optional
        TCR chain: alpha, beta or paired. The default is 'beta'.
    sparse : Bool, optional
        Turn on sparse distance computation. The default is False.

    Returns
    -------
    S : numpy.array
        TCRDist distance matrix.
    seq : pandas.Series
        pandas.Series with sequences for which distances have been calculated.
    gt : pandas.DataFrame
        Ground truth. pandas.DataFrame containing information about the
        TCR sequence and its cognate epitope target.

    '''   
    if data is None:    
        vdjdb = parse_vdjdb(os.path.abspath('./clustcr/input/vdjdb/vdjdb_full.txt'), q=1)
    else:
        vdjdb = data

    if chain == 'beta':
        cdr3 = 'cdr3_b_aa'
        v_name = 'v_b_gene'
        vdjdb = vdjdb.drop(columns=['cdr3.alpha', 'v.alpha'])
        vdjdb = vdjdb.rename(columns={'cdr3.beta':cdr3,
                                      'v.beta':v_name})
    elif chain == 'alpha':
        cdr3 = 'cdr3_a_aa'
        v_name = 'v_a_gene'
        vdjdb = vdjdb.drop(columns=['cdr3.beta', 'v.beta'])
        vdjdb = vdjdb.rename(columns={'cdr3.alpha':cdr3,
                                      'v.alpha':v_name})        
    
    df_epi = vdjdb[[cdr3, v_name, 'antigen.epitope']].dropna().drop_duplicates()
    seq = df_epi.drop(columns = ['antigen.epitope']).drop_duplicates().reset_index(drop=True)
    
    gt = df_epi.rename(columns = {cdr3:'CDR3',
                                  v_name:'V',
                                  'antigen.epitope':'Epitope'})
    
    if sparse:
        
        tr = TCRrep(cell_df = seq,
                    organism = 'human',
                    chains = ['beta'],
                    db_file = 'alphabeta_gammadelta_db.tsv',
                    compute_distances = False)
        
        tr.cpus = 2
        tr.compute_sparse_rect_distances(radius = 200, chunk_size = 500)
        S = tr.rw_beta
        
    else:
        
        tr = TCRrep(cell_df = seq,
            organism = 'human',
            chains = [chain],
            db_file = 'alphabeta_gammadelta_db.tsv',
            compute_distances = True)
    
        S = tr.pw_cdr3_b_aa
        
    return S, seq, gt

def normalize(edge):
    '''
    Introduce normalization property on an edge that is represented as a tuple.
    The normalization property will constraint the ordering of two nodes that
    make up an edge. This prevents duplicated edges.

    Parameters
    ----------
    edge : tuple
        Tuple of length two, containing two nodes, represented as integers.

    Returns
    -------
    (n1, n2) : tuple
        Normalized edge.
        
    '''
    n1, n2 = edge
    if n1 > n2:
        n1, n2 = n2, n1
    return (n1, n2)

def greedy_clustering(dm, threshold):
    '''
    Greedy graph clustering approach that uses a fixed distance-threshold to 
    assign nodes to cluster. The algorithm starts by computing all pairs 
    of sequences that satisfy the predefined distance threshold (edge list). 
    Next, it finds the sequence with the highest degree (i.e. the most neighbors), 
    assigns this as the cluster centre, and removes it and its neighbors 
    from the edge list. This process is repeated until all sequences are clustered.

    Parameters
    ----------
    dm : numpy.array
        Distance matrix.
    threshold : int
        Distance threshold for defining network edges.

    Returns
    -------
    res : pandas.DataFrame
        Dataframe containing clustering results.

    '''

    edges = np.argwhere(dm<=threshold)
    print(len(edges))
    edges = set(map(normalize, edges)) # Remove duplicated edges
    edges = np.array(list(edges)) # Edgelist to array
    print(len(edges))
    
    cid = 0
    res = pd.DataFrame()
    
    while len(edges) > 0:
        
        G = nx.from_edgelist(edges)
        degrees = pd.DataFrame(G.degree(), columns=['node', 'degree'])
        degrees = degrees.set_index('node')
        degrees = degrees.sort_values(by='degree', ascending=False)
        max_degree = degrees.idxmax().values
        
        cluster = edges[np.where(edges[:,0]==max_degree)]
        ids = np.unique(cluster)
        cids = [cid] * len(ids)

        if len(ids) <= 1:
            break
        
        cluster_iter = pd.DataFrame({'seq_id':ids,'cluster':cids})
        res = res.append(cluster_iter)
        
        for i in ids:
            edges = edges[np.where(edges[:,0]!=i)] # Remove from column 1
            edges = edges[np.where(edges[:,1]!=i)] # Remove from column 2
        
        cid += 1
            
    return res

def cluster_TCRDist_matrix(dm=None, cdr3=None, gt=None, method='DBSCAN'):
    '''
    Function for clustering of the TCRDist-calculated distance matrix.
    The function provides several methods for clustering, which include:
        - Greedy: fixed-distance threshold clustering method that groups
        of sequences that are connected in a graph.
        - DBSCAN: density-based clustering method that groups of densely
        packed points.
        - Kmeans: iterative clustering approach that partitions the data
        into k clusters.

    Parameters
    ----------
    dm : numpy.array, optional
        TCRDist-calculated distance matrix. The default is None.
    cdr3 : pandas.Series, optional
        pandas.Series containing the input CDR3 sequences. The default is None.
    gt : pandas.DataFrame, optional
        Ground truth data. The default is None.
    method : str, optional
        Clustering method used to cluster the TCRDist distance matrix. 
        Available methods include: greedy, DBSCAN and Kmeans.
        The default is 'DBSCAN'.

    Returns
    -------
    Clustering results

    '''
    # Available methods
    methods = ['GREEDY', 'DBSCAN', 'KMEANS']
    assert method.upper() in methods, r'Please choose one of the following: /n %s' % methods
    
    # If any of the parameters not specified, compute it using default settings
    if dm is None:
        dm, cdr3, gt = TCRDist(sparse=False)
    if gt is None:
        dm, cdr3, gt = TCRDist(sparse=False)
    if cdr3 is None:
        dm, cdr3, gt = TCRDist(sparse=False)
             
    if method.upper() == 'GREEDY':    
        # Greedy clustering
        clusters = greedy_clustering(dm, 12)
        clusters = clusters.rename(columns={'seq_id':'Index'})
        clusters = clusters.set_index('Index', drop=True)
        clusters = clusters.merge(right=cdr3, left_index=True, right_index=True)
        clusters = clusters.rename(columns={'cdr3_b_aa':'CDR3',
                                            'v_b_gene':'V'})
        metrics = ClusteringResult(clusters).metrics(gt)
        return metrics.summary()

    elif method.upper() == 'DBSCAN':
        # DBSCAN
        clustering = DBSCAN(eps=100, min_samples=2, n_jobs=-1).fit(dm)
        labels = clustering.labels_
        clusters = cdr3.rename(columns={'cdr3_b_aa':'CDR3',
                                        'v_b_gene':'V'})
        clusters['cluster'] = labels
        clusters = clusters[clusters['cluster']!=-1]
        metrics = ClusteringResult(clusters).metrics(gt)
        return metrics.summary()
    
    else:
        # K-Means
        kmeans = KMeans(n_clusters=500).fit(dm)
        labels = kmeans.labels_
        clusters = cdr3.rename(columns={'cdr3_b_aa':'CDR3',
                                        'v_b_gene':'V'})
        clusters['cluster'] = labels
        metrics = ClusteringResult(clusters).metrics(gt)
        return metrics.summary()