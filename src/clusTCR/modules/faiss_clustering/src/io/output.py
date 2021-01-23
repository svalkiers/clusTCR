import pandas as pd
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, 'results')


def path_in_results(filename):
    return os.path.join(RESULTS, filename)


def simple_cluster_output(clustering, data):
    s = ''
    for id_list in clustering.get_cluster_contents():
        s += simple_single_cluster_output(id_list, data)
    return s


def simple_single_cluster_output(cluster, data: pd.Series):
    """
    Takes a list of IDs and the data
    Returns a simple output of the clustered sequences
    """
    s = '-----------------------------------------\n'
    for id in cluster:
        s += str(data[id]) + '\n'
    return s


def simple_output_to_file(index, data: pd.Series, fname):
    if fname is None:
        return
    with open(fname, 'w') as f:
        f.write(simple_cluster_output(index, data))


def cluster_to_csv(clustering, data, fname):
    cluster_contents_to_csv(clustering.get_cluster_contents(), data, fname)


def cluster_contents_to_csv(contents, data, fname):
    s = ''
    for cluster in contents:
        for id in cluster:
            s += str(data[id])
            if id != cluster[-1]:
                s += ','
        s += '\n'
    with open(fname, 'w') as f:
        f.write(s)

