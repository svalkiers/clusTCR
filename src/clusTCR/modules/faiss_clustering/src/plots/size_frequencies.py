from ..cluster import *
from ..io.output import get_cluster_contents
import matplotlib.pyplot as plt


def main():
    cdr3, epitopes = small_vdj()
    size_frequencies(cdr3)


def size_frequencies(data):
    for combo in PROPERTY_COMBINATIONS:
        index = profile_cluster(data, combo)
        name = 'cluster by ' + ' + '.join(combo)
        size_frequencies_plot(index, name)


def size_frequencies_plot(index, name):
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
    plt.show()


if __name__ == '__main__':
    main()
