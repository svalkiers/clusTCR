import pandas as pd
import networkx as nx
import markov_clustering as mcl
import parmap
import multiprocessing
from Bio import pairwise2
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio.SubsMat import MatrixInfo
from clustcr.clustering.tools import create_edgelist


def del_list_inplace(l, id_to_del):
    for i in sorted(id_to_del, reverse=True):
        del(l[i])


def blosum_cutoff(edgelist, cutoff):
    edgelist = list(edgelist)
    to_delete = []
    matrix = MatrixInfo.blosum62
    for i in range(len(edgelist)):
        a, b = edgelist[i].split()
        max1 = pairwise2.align.globaldx(a, a, matrix, score_only=True)
        max2 = pairwise2.align.globaldx(b, b, matrix, score_only=True)
        maximum = max(max1, max2)
        score = pairwise2.align.globaldx(a, b, matrix, score_only=True)
        score /= maximum
        if score < cutoff:
            to_delete.append(i)
    del_list_inplace(edgelist, to_delete)
    return edgelist


def MCL(cdr3, edgelist=None, bl_cutoff=0.93, mcl_hyper=[1.2, 2], outfile=None):
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

    clusters = {"CDR3": [], "cluster": []}

    # G = nx.parse_adjlist(edgelist, nodetype=str)
    if bl_cutoff:
        edgelist = blosum_cutoff(edgelist, bl_cutoff)
    if not edgelist:
        return pd.DataFrame(data=clusters)
    G = nx.Graph()
    for elem in edgelist:
        a, b = elem.split()
        G.add_edge(a, b)
    m = nx.to_scipy_sparse_matrix(G)

    # Run MC
    result = mcl.run_mcl(m, inflation=mcl_hyper[0], expansion=mcl_hyper[1])
    mcl_output = mcl.get_clusters(result)
    identifiers = list(G.nodes())

    # Map cluster ids back to seqs
    cluster_ids = dict()
    count = 0
    for i in range(len(mcl_output)):
        if len(mcl_output[i]) <= 1:
            continue
        cluster_ids[count] = list(identifiers[node] for node in mcl_output[i])
        count += 1

    # Generate nodelist
    for c in cluster_ids:
        for seq in cluster_ids[c]:
            clusters["CDR3"].append(seq)
            clusters["cluster"].append(c)
    clusters = pd.DataFrame(data=clusters)

    # Write to file
    if outfile is not None:
        clusters.to_csv(outfile, sep="\t", index=False)

    return clusters


def MCL_multi(edgelist, cdr3, bl_cutoff):
    return MCL(cdr3, edgelist, bl_cutoff)


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


def MCL_multiprocessing_from_preclusters(cdr3, preclust, n_cpus, blosum_cutoff, hd_cutoff):
    """
    Pool multiple processes for parallelization using multiple cpus.
    """
    import time
    a = time.time()
    cluster_contents = preclust.cluster_contents()
    edges = {i: create_edgelist(cluster, hd_cutoff) for i, cluster in enumerate(cluster_contents)}
    remaining_edges = edges.values()
    print(time.time() - a, 'edgelist')
    # Perform MCL on other clusters
    with multiprocessing.Pool(n_cpus) as pool:
        nodelist = parmap.map(MCL_multi,
                              remaining_edges,
                              cdr3,
                              blosum_cutoff,
                              pm_parallel=True,
                              pm_pool=pool)
    # Fix cluster ids
    for c in range(len(nodelist)):
        if c != 0:
            nodelist[c]['cluster'] += nodelist[c - 1]['cluster'].max() + 1
    return pd.concat(nodelist, ignore_index=True)


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
