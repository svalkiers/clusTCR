from src.cluster import *
from src.io.datasets import covid19_repertoire
from src.io.output import get_cluster_contents, simple_single_cluster_output, path_in_results
import matplotlib.pyplot as plt


def main():
    data = covid19_repertoire()
    size_frequencies(data)
    random_representative_comparison(data)


def random_representative_comparison(data):
    representative = data.sample().keys()[0]
    s = ''
    for combo in PROPERTY_COMBINATIONS:
        index = profile_cluster(data, combo)
        cluster = get_cluster_of_representative(representative, index)
        s += ' + '.join(combo).upper() + '\n'
        s += output_single_cluster_with_representative(cluster, representative, data) + '\n'
    with open(path_in_results('profile_method/comparison/representative_clusters.txt'), 'w') as f:
        f.write(s)


def output_single_cluster_with_representative(cluster, representative, data):
    output = simple_single_cluster_output(cluster, data)
    representative_sequence = str(data[representative])
    return output.replace(representative_sequence, '> ' + representative_sequence)


def get_cluster_of_representative(representative, index):
    for cluster in get_cluster_contents(index):
        if representative in cluster:
            return cluster
    return None


def size_frequencies(data):
    for combo in PROPERTY_COMBINATIONS:
        index = profile_cluster(data, combo)
        fname = '_'.join(combo)
        name = 'cluster by ' + ' + '.join(combo)
        size_frequencies_plot(index, name, fname)


def size_frequencies_plot(index, name, fname):
    contents = get_cluster_contents(index)
    sizes = list(map(lambda x: len(x), contents))
    frequencies = {}
    for size in sizes:
        frequencies[size] = frequencies.get(size, 0) + 1
    x = sorted(frequencies.keys())
    y = []
    for point in x:
        y.append(frequencies[point])
    plt.figure()
    plt.bar(x, y)
    plt.axis([0, 50, 0, 300])
    plt.title(name)
    plt.xlabel('cluster size')
    plt.ylabel('amount of clusters')
    plt.savefig(fname=path_in_results(f'profile_method/comparison/{fname}.png'))


if __name__ == '__main__':
    main()
