import pandas as pd
import networkx as nx
import markov_clustering as mcl
import parmap
import multiprocessing
from clustcr.clustering.tools import create_edgelist


def MCL(cdr3, edgelist=None, mcl_hyper=[1.2, 2], outfile=None):
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
        edgelist = create_edgelist(cdr3)

    G = nx.parse_adjlist(edgelist, nodetype=str)
    m = nx.to_scipy_sparse_matrix(G)

    # Run MC
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

    return clusters


def MCL_multi(edgelist, cdr3):
    return MCL(cdr3, edgelist)


def MCL_multiprocessing_from_preclusters(cdr3, preclust, n_cpus):
    # Pool multiple processes for parallelization using multiple cpus.
    clusters = []
    idxs_to_remove = []
    cluster_contents = preclust.cluster_contents()
    edges = [create_edgelist(c) for c in cluster_contents]
    # Clusters containing no edges with HD = 1 are isolated
    for val in edges:
        if len(val) == 0:
            idx = edges.index(val)
            idxs_to_remove.append(idx)
            clust = cluster_contents[idx]
            if len(clusters) == 0:
                clusters.append(pd.DataFrame({'CDR3': clust,
                                              'cluster': [0] * len(clust)}))
            else:
                clusters.append(pd.DataFrame({'CDR3': clust,
                                              'cluster': (clusters[-1]['cluster'].max() + 1) * len(clust)}))
    for index in sorted(idxs_to_remove, reverse=True):
        del edges[index]
    # Perform MCL on other clusters
    with multiprocessing.Pool(n_cpus) as pool:
        nodelist = parmap.map(MCL_multi,
                              edges,
                              cdr3,
                              pm_parallel=True,
                              pm_pool=pool)
        nodelist += clusters
    # Fix cluster ids
    for c in range(len(nodelist)):
        if c == 0:
            pass
        else:
            nodelist[c]['cluster'] += nodelist[c - 1]['cluster'].max() + 1
    return pd.concat(nodelist)


def MCL_from_preclusters(cdr3, preclust):
    initiate = True
    for c in preclust.cluster_contents():
        try:
            edges = create_edgelist(c)
            if initiate:
                nodelist = MCL(cdr3, edges)
                initiate = False
            else:
                nodes = MCL(cdr3, edges)
                nodes["cluster"] = nodes["cluster"] + nodelist["cluster"].max() + 1
                nodelist = nodelist.append(nodes)
        # If no edges can be found, leave cluster as is
        except nx.NetworkXError:
            cluster = pd.DataFrame({"CDR3": c,
                                    "cluster": [nodelist["cluster"].max() + 1] * len(c)})
            nodelist = nodelist.append(cluster)
    return nodelist
