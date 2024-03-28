import pandas as pd
import networkx as nx
import markov_clustering as mcl
import community
import parmap
import multiprocessing
from clustcr.clustering.tools import create_edgelist

def MCL(cdr3=None, edgelist=None, mcl_hyper=[1.2, 2], outfile=None):
    """
    Perform clustering on a network of CDR3 amino acid sequences with
    a known hamming distance, using the Markov clustering (MCL) algorithm.
    For more info about the inflation and expansion parameters,
    visit: https://micans.org/mcl/


    Parameters
    ----------
    edgelist : set, optional
        Tab-separated edgelist. The default is None.
    mcl_hyper : list, optional
        MCL hyperparameters: inflation and expansion.
        The default is [1.2,2].
    outfile : str, optional
        Name of outfile. The default is None.

    Returns
    -------
    clusters : pd.DataFrame
        pd.DataFrame containing two columns: 'junction_aa' and 'cluster'.
        The first column contains CDR3 sequences, the second column
        contains the corresponding cluster ids.
    """
    if edgelist is None:
        assert cdr3 is not None, 'No sequences or edges provided.'
        edgelist = create_edgelist(cdr3)

    try:
        G = nx.parse_adjlist(edgelist, nodetype=str)
        m = nx.to_scipy_sparse_array(G)
    
        # Run MCL
        result = mcl.run_mcl(m, inflation=mcl_hyper[0], expansion=mcl_hyper[1])
        mcl_output = mcl.get_clusters(result)
        identifiers = list(G.nodes())
    
        # Map cluster ids back to seqs
        cluster_ids = dict()
        for i in range(len(mcl_output)):
            cluster_ids[i] = list(identifiers[i] for i in mcl_output[i])
    
        # Generate nodelist
        clusters = {"junction_aa": [], "cluster": []}
        for c in cluster_ids:
            for seq in cluster_ids[c]:
                clusters["junction_aa"].append(seq)
                clusters["cluster"].append(c)
        clusters = pd.DataFrame(data=clusters)
    
        # Write to file
        if outfile is not None:
            clusters.to_csv(outfile, sep="\t", index=False)
    except nx.NetworkXError:
        clusters = pd.DataFrame({"junction_aa": [], "cluster": []})

    return clusters


def louvain(cdr3, edgelist=None):
    if edgelist is None:
        edgelist = create_edgelist(cdr3)

    try:
        G = nx.parse_adjlist(edgelist, nodetype=str)
        partition = community.best_partition(G)
    except nx.NetworkXError:
        partition = pd.DataFrame({"junction_aa": [], "cluster": []})
        
    return pd.DataFrame.from_dict(
        partition, orient="index", columns=["cluster"]
        ).reset_index().rename(columns={'index': 'junction_aa'})

def MCL_multi(edgelist=None, cdr3=None, mcl_hyper=[1.2,2]):
    return MCL(cdr3, edgelist, mcl_hyper=mcl_hyper)


def clusters_without_hd1_edges(edges, cluster_contents):
    """
    Returns clusters that don't contain edges with edit distance 1
    Also removes them from the edges
    """
    clusters = []
    ids_to_be_removed = []
    for i, edge_list in edges.items():
        if len(edge_list) != 0:
            continue
        ids_to_be_removed.append(i)
        cluster = cluster_contents[i]
        clusters.append(pd.DataFrame({
            'junction_aa': cluster,
            'cluster': [0] * len(cluster)
        }))
    for id in ids_to_be_removed:
        del edges[id]
    return clusters

def clean_edgelist(edges):
    ids_to_be_removed = []
    for i, edge_list in edges.items():
        if len(edge_list) != 0:
            continue
        ids_to_be_removed.append(i)
    clean_edges = {i:edges[i] for i in edges.keys() if i not in ids_to_be_removed}
    return clean_edges

def MCL_multiprocessing_from_preclusters(preclust, mcl_hyper, n_cpus):
    """
    Pool multiple processes for parallelization using multiple cpus.
    """
    cluster_contents = preclust.cluster_contents()
    edges = {i: create_edgelist(cluster) for i, cluster in enumerate(cluster_contents)}
    # Clusters containing no edges with HD = 1 are removed
    clean_edges = clean_edgelist(edges)
    remaining_edges = clean_edges.values()
    # Perform MCL on other clusters
    cdr3 = None
    with multiprocessing.Pool(n_cpus) as pool:
        nodelist = parmap.map(MCL_multi,
                              remaining_edges,
                              cdr3,
                              mcl_hyper=mcl_hyper,
                              pm_parallel=True,
                              pm_pool=pool)

    # Fix cluster ids
    for c in range(len(nodelist)):
        if c != 0:
            nodelist[c]['cluster'] += nodelist[c - 1]['cluster'].max() + 1
    return pd.concat(nodelist, ignore_index=True)

def MCL_multiprocessing_from_preclusters_test(preclust, mcl_hyper, n_cpus):
    """
    BUGFIXING, DON'T USE
    """
    cluster_contents = preclust.cluster_contents()
    edges = {i: create_edgelist(cluster) for i, cluster in enumerate(cluster_contents)}
    # Clusters containing no edges with HD = 1 are isolated
    # clusters = clusters_without_hd1_edges(edges, cluster_contents)
    remaining_edges = edges.values()
    cdr3 = None
    # Perform MCL on other clusters
    with multiprocessing.Pool(n_cpus) as pool:
        nodelist = parmap.map(MCL_multi,
                              remaining_edges,
                              cdr3,
                              mcl_hyper=mcl_hyper,
                              pm_parallel=True,
                              pm_pool=pool)
        # nodelist += clusters

    # Fix cluster ids
    for c in range(len(nodelist)):
        if c != 0:
            nodelist[c]['cluster'] += nodelist[c - 1]['cluster'].max() + 1
    return pd.concat(nodelist, ignore_index=True)

def MCL_from_preclusters(preclust, mcl_hyper):
    '''
    Second pass using MCL, without multiprocessing.
    '''
    initiate = True
    nodelist = pd.DataFrame(columns=["junction_aa", "cluster"])
    for c in preclust.cluster_contents():
        edges = create_edgelist(c)
        if initiate:
            nodes = MCL(edgelist=edges, mcl_hyper=mcl_hyper)
            nodelist = pd.concat([nodelist,nodes],ignore_index=True)
            initiate = False
        else:
            nodes = MCL(edgelist=edges, mcl_hyper=mcl_hyper)
            nodes["cluster"] = nodes["cluster"] + nodelist["cluster"].max() + 1
            nodelist = pd.concat([nodelist,nodes],ignore_index=True)
    return nodelist

def louvain_multiprocessing_from_preclusters(cdr3, preclust, n_cpus):
    """
    Pool multiple processes for parallelization using multiple cpus.
    """
    cluster_contents = preclust.cluster_contents()
    edges = {i: create_edgelist(cluster) for i, cluster in enumerate(cluster_contents)}
    # Clusters containing no edges with HD = 1 are isolated
    clean_edges = clean_edgelist(edges)
    remaining_edges = clean_edges.values()
    # Perform MCL on other clusters
    with multiprocessing.Pool(n_cpus) as pool:
        nodelist = parmap.map(louvain,
                              remaining_edges,
                              cdr3,
                              pm_parallel=True,
                              pm_pool=pool)
    # Fix cluster ids
    for c in range(len(nodelist)):
        if c != 0:
            nodelist[c]['cluster'] += nodelist[c - 1]['cluster'].max() + 1
    return pd.concat(nodelist, ignore_index=True)

def louvain_from_preclusters(cdr3, preclust):
    '''
    Second pass using Louvain clustering, without multiprocessing.
    '''
    initiate = True
    nodelist = pd.DataFrame(columns=["junction_aa", "cluster"])
    for c in preclust.cluster_contents():
        edges = create_edgelist(c)
        if initiate:
            nodes = louvain(cdr3, edges)
            nodelist = pd.concat([nodelist,nodes],ignore_index=True)
            initiate = False
        else:
            nodes = louvain(cdr3, edges)
            nodes["cluster"] = nodes["cluster"] + nodelist["cluster"].max() + 1
            nodelist = pd.concat([nodelist,nodes],ignore_index=True)
    return nodelist
