import json
import pandas as pd

from os.path import dirname, abspath, join

ROOT = dirname(dirname(dirname(abspath(__file__))))
DATA = join(ROOT, 'data')

def path_in_data(filename):
    return join(DATA, filename)

def from_tcrdata(filename, epitopes = False, sep = '\t', q = None):
    # NOTE: q-score is only used for VDJdb data.
    bm = pd.read_csv(filename, sep = sep)
    if q is not None:
        bm = bm[bm["Score"] >= q]
    if epitopes:
        return bm[["CDR3", "Epitope"]].reset_index(drop = True)
    else:
        return bm["CDR3"].unique()
    
def from_edgelist(filename):
    
    with open(filename, 'r') as f:
        edges = set(f.read().splitlines())
    
    return edges
    
def from_json(filename, max_distance = 3, weight = -3):
    '''
    Import your weighted network.
    
    jsonfile : jsonfile that contains the 
    max_distance : max allowed distance
    weight : negative exponential that downweighs sequence pairs proportional to their distance
    '''
    
    with open(filename, 'r') as f:
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
            edges.append((*tuple(dist[key][val]), int(key)**(weight)))
        
    return edges