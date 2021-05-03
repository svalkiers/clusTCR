import pandas as pd
import networkx as nx
import markov_clustering as mcl
import parmap
import multiprocessing
import time
from clustcr.clustering.tools import create_edgelist

def timeit(myfunc):
    # Decorator to keep track of time required to run a function
    def timed(*args, **kwargs):
        start = time.time()
        result = myfunc(*args, **kwargs)
        end = time.time()
        print(f'Total time to run \'{myfunc.__name__}\': {(end-start):.3f}s')
        return result
    return timed


def MCL(cdr3, edgelist=None, distance_metric='HAMMING', mcl_hyper=[1.2, 2], outfile=None):
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
        pd.DataFrame containing two columns: 'CDR3' and 'cluster'.
        The first column contains CDR3 sequences, the second column
        contains the corresponding cluster ids.
    """
    if edgelist is None:
        edgelist = create_edgelist(cdr3, method=distance_metric)

    try:
        G = nx.parse_adjlist(edgelist, nodetype=str)
        m = nx.to_scipy_sparse_matrix(G)
    
        # Run MCL
        result = mcl.run_mcl(m, inflation=mcl_hyper[0], expansion=mcl_hyper[1])
        mcl_output = mcl.get_clusters(result)
        identifiers = list(G.nodes())
    
        # Map cluster ids back to seqs
        cluster_ids = dict()
        for i in range(len(mcl_output)):
            cluster_ids[i] = list(identifiers[i] for i in mcl_output[i])
    
        # Generate nodelist
        clusters = {"CDR3": [], "cluster": []}
        for c in cluster_ids:
            for seq in cluster_ids[c]:
                clusters["CDR3"].append(seq)
                clusters["cluster"].append(c)
        clusters = pd.DataFrame(data=clusters)
    
        # Write to file
        if outfile is not None:
            clusters.to_csv(outfile, sep="\t", index=False)
    except nx.NetworkXError:
        clusters = pd.DataFrame({"CDR3": [], "cluster": []})

    return clusters


def MCL_multi(edgelist, cdr3, mcl_hyper=[1.2,2]):
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
            'CDR3': cluster,
            'cluster': [0] * len(cluster)
        }))
    for id in ids_to_be_removed:
        del edges[id]
    return clusters

def MCL_multiprocessing_from_preclusters(cdr3, preclust, distance_metric, mcl_hyper, n_cpus):
    """
    Pool multiple processes for parallelization using multiple cpus.
    """
    cluster_contents = preclust.cluster_contents()
    edges = {i: create_edgelist(cluster, method=distance_metric) for i, cluster in enumerate(cluster_contents)}
    # Clusters containing no edges with HD = 1 are isolated
    clusters = clusters_without_hd1_edges(edges, cluster_contents)
    remaining_edges = edges.values()
    # Perform MCL on other clusters
    with multiprocessing.Pool(n_cpus) as pool:
        nodelist = parmap.map(MCL_multi,
                              remaining_edges,
                              cdr3,
                              mcl_hyper=mcl_hyper,
                              pm_parallel=True,
                              pm_pool=pool)
        nodelist += clusters

    # Fix cluster ids
    for c in range(len(nodelist)):
        if c != 0:
            nodelist[c]['cluster'] += nodelist[c - 1]['cluster'].max() + 1
    return pd.concat(nodelist, ignore_index=True)

def MCL_from_preclusters(cdr3, preclust, distance_metric, mcl_hyper):
    initiate = True
    nodelist = pd.DataFrame(columns=["CDR3", "cluster"])
    for c in preclust.cluster_contents():
        try:
            edges = create_edgelist(c, method=distance_metric)
            if initiate:
                nodes = MCL(cdr3, edges, mcl_hyper=mcl_hyper)
                nodelist = nodelist.append(nodes)
                initiate = False
            else:
                nodes = MCL(cdr3, edges, mcl_hyper=mcl_hyper)
                nodes["cluster"] = nodes["cluster"] + nodelist["cluster"].max() + 1
                nodelist = nodelist.append(nodes)
        # If no edges can be found, leave cluster as is
        except nx.NetworkXError:
            try:
                cluster = pd.DataFrame({"CDR3": c,
                                        "cluster": [nodelist["cluster"].max() + 1] * len(c)})
                nodelist = nodelist.append(cluster)
            except KeyError:
                cluster = pd.DataFrame({"CDR3": c, "cluster": [0] * len(c)})
                nodelist = nodelist.append(cluster)
    return nodelist
