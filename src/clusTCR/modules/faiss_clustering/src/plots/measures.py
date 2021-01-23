from ..evaluate import calculate_consistency, calculate_purity, calculate_single_cluster_purity
from ..cluster import *
from ..io.datasets import *
from ..io.output import get_cluster_contents
import matplotlib.pyplot as plt


def main():
    cdr3, epitopes = small_vdj()
    plot_measure_with_prop_combinations(cdr3, epitopes,
                                        measure_func=calculate_purity,
                                        title='VDJdb subset')
    plot_measure_with_prop_combinations(cdr3, epitopes,
                                        measure_func=calculate_consistency,
                                        title='VDJdb subset')


def bar_chart_purity(cdr3, epitopes, title, remove_singular_clusters=True):
    for i in [3, 5, 7, 9, 11, 13]:
        plt.figure()

        index = profile_cluster(cdr3, ALL_PROPS, items_per_cluster=i)
        counts = {}

        for cluster in get_cluster_contents(index):
            if remove_singular_clusters and len(cluster) <= 1:
                continue
            cluster_purity = calculate_single_cluster_purity(cluster, epitopes)
            percentage = int(cluster_purity * 100)
            counts[percentage] = counts.get(percentage, 0) + 1

        x = counts.keys()
        y = counts.values()

        plt.bar(x, y)
        plt.xlabel('Purity')
        plt.ylabel('Amount of clusters')
        plt.title(f'{title} (avg {i} TCRs per cluster)')
        plt.show()


def plot_purity_vs_consistency(cdr3, epitopes, title):
    fig, (ax1, ax2) = plt.subplots(2)

    x = []
    y_consistency = []
    y_purity = []

    for i in [3, 5, 7, 9, 11, 13]:
        index = profile_cluster(cdr3, ALL_PROPS, items_per_cluster=i)
        x.append(i)
        y_purity.append(calculate_purity(index, epitopes, cdr3.size))
        y_consistency.append(calculate_consistency(index, epitopes, cdr3.size))

    ax1.plot(x, y_consistency, label='consistency')
    ax1.set(ylabel='consistency')

    ax2.plot(x, y_purity, label='purity')
    ax2.set(ylabel='purity')

    fig.suptitle(title)
    plt.xlabel("Average amount of TCR's per cluster")
    plt.show()


def plot_measure_with_prop_combinations(cdr3, epitopes, measure_func, title):
    plt.figure()

    COMBOS = [
        [BASICITY, HELICITY, HYDROPHOBICITY],
        Z_SCORES,
        [ISOELECTRIC, *Z_SCORES],
        [MUTATION_STABILITY, *Z_SCORES],
        [BASICITY, HELICITY, *Z_SCORES],
        [BASICITY, HELICITY, HYDROPHOBICITY, *Z_SCORES],
        [BASICITY, HELICITY, HYDROPHOBICITY, MUTATION_STABILITY],
        [BASICITY, HELICITY, MUTATION_STABILITY, *Z_SCORES],
        [BASICITY, HELICITY, MUTATION_STABILITY, ISOELECTRIC, *Z_SCORES],
    ]

    for combo in COMBOS:
        x = []
        y = []
        for i in range(1, 20):
            index = profile_cluster(cdr3, combo, items_per_cluster=i)
            x.append(i)
            y.append(measure_func(index, epitopes, cdr3.size))
        atchley_string = ' + '.join(ATCHLEY_FACTORS)
        z_string = ' + '.join(Z_SCORES)
        plt.plot(x, y, label=' + '.join(combo).replace(atchley_string, 'atchley').replace(z_string, 'z_scores'))

    plt.xlabel("Average amount of TCR's per cluster")
    plt.ylabel(measure_func.__name__.title())

    plt.title(title)
    plt.legend(loc='best')
    plt.show()


def plot_measure_with_prop_combinations_and_compare(cdr3, epitopes, measure_func, title):
    plt.figure()

    COMBOS = [
        [BASICITY, HELICITY, HYDROPHOBICITY],
        [ISOELECTRIC],
        [MUTATION_STABILITY],
        [BASICITY, HELICITY, HYDROPHOBICITY, ISOELECTRIC],
        [BASICITY, HELICITY, HYDROPHOBICITY, MUTATION_STABILITY],
        [BASICITY, HELICITY, HYDROPHOBICITY, MUTATION_STABILITY, ISOELECTRIC],
    ]

    ax = plt.gca()
    for combo in COMBOS:
        x = []
        y1 = []
        y2 = []
        for i in range(1, 20):
            index1 = profile_cluster(cdr3, combo, items_per_cluster=i, add_average=False)
            index2 = profile_cluster(cdr3, combo, items_per_cluster=i, add_average=True)
            x.append(i)
            y1.append(measure_func(index1, epitopes, cdr3.size))
            y2.append(measure_func(index2, epitopes, cdr3.size))
        color = next(ax._get_lines.prop_cycler)['color']
        plt.plot(x, y2, label=' + '.join(combo), color=color)
        plt.plot(x, y1, linestyle=':', color=color)

    plt.xlabel("Average amount of TCR's per cluster")
    plt.ylabel(measure_func.__name__.title())

    plt.title(title)
    plt.legend(loc='best')
    plt.show()


if __name__ == '__main__':
    main()
