import json
import pandas as pd

from os.path import dirname, abspath, join

DIR = dirname(abspath(__file__))
ROOT = dirname(dirname(dirname(abspath(__file__))))
DATA = join(ROOT, 'data')


def path_in_data(filename):
    return join(DATA, filename)


def imgt_v_genes(filename='alphabeta_gammadelta_db.tsv'):
    v_genes = pd.read_csv(join(DIR, filename), sep='\t')
    v_genes = v_genes[v_genes['chain'] == 'B']
    v_genes = v_genes[v_genes['region'] == 'V']
    v_genes = v_genes[v_genes['organism'] == 'human']
    return v_genes['id']


def from_edgelist(filename):
    with open(filename, 'r') as f:
        edges = set(f.read().splitlines())

    return edges


def weighted_network_from_json(filename, max_distance=3, weight=-1):
    """
    Import weighted network from json file.
    max_distance : max allowed distance
    weight : negative exponential that downweighs sequence pairs proportional to their distance
    """

    with open(path_in_data(filename), 'r') as f:
        dist = json.loads(f.read())

    distances = []
    for d in dist.keys():
        distances.append(int(d))

    for i in range(max_distance + 1, max(distances) + 1, 1):
        try:
            del dist[str(i)]
        except KeyError:
            pass

    edges = []
    for key in dist.keys():
        for val in range(len(dist[key])):
            edges.append((*tuple(dist[key][val]), int(key) ** (weight)))

    return edges
