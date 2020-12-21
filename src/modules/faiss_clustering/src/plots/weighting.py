from ..evaluate import calculate_purity
from ..cluster import *
from ..io.datasets import vdj
import matplotlib.pyplot as plt


def main():
    cdr3, epitopes = vdj()
    plot_weighted_measure(cdr3, epitopes,
                          title='VDJdb',
                          measure_func=calculate_purity,
                          weighting_func=weight_vector_center,
                          weighting_name='Center')


def weight_vector_sides(vec, factor_percentage, percentage_of_vec_weighed):
    """ Factor and amount between 0 and 1 """
    n = len(vec)
    n_transform = int(n * percentage_of_vec_weighed / 2)
    for i in range(n_transform):
        vec[i] = vec[i] * factor_percentage
        vec[- i - 1] = vec[- i - 1] * factor_percentage
    return vec


def weight_vector_center(vec, factor_percentage, percentage_of_vec_weighed):
    """ Factor and amount between 0 and 1 """
    n = len(vec)
    n_transform = int(n * percentage_of_vec_weighed / 2)
    middle = n // 2
    for i in range(n_transform):
        left = middle - 1 - i
        right = middle + i
        vec[left] = vec[left] * factor_percentage
        vec[right] = vec[right] * factor_percentage
    return vec


def test_weighted_clustering(cdr3, epitopes):
    items_per_cluster = 7
    index1 = profile_cluster(cdr3, items_per_cluster=items_per_cluster)
    purity1 = calculate_purity(index1, epitopes, cdr3.size)
    index2 = profile_cluster(cdr3, items_per_cluster=items_per_cluster,
                             vector_mapping_func=lambda vec: weight_vector_center(vec, 0.4, 0.15))
    purity2 = calculate_purity(index2, epitopes, cdr3.size)

    print(purity1)
    print(purity2)


def plot_weighted_measure(cdr3, epitopes, title, measure_func, weighting_func, weighting_name, items_per_cluster=7):
    plt.figure()
    for factor in np.arange(0.1, 1, 0.1):
        x = []
        y = []
        for vec_percentage in np.arange(0, 1, 0.05):
            vec_mapping = lambda vec: weighting_func(vec, factor, vec_percentage)
            index = profile_cluster(cdr3, properties=[BASICITY, HELICITY, HYDROPHOBICITY],
                                    items_per_cluster=items_per_cluster, vector_mapping_func=vec_mapping)
            index_purity = measure_func(index, epitopes, cdr3.size)
            x.append(vec_percentage)
            y.append(index_purity)
        plt.plot(x, y, label=f'*{factor:.2f}')

    plt.xlabel('Percentage of vector weighted')
    plt.ylabel(measure_func.__name__.title())
    plt.legend(loc='best')
    plt.title(f'{weighting_name} weighted of {title}')
    plt.show()


if __name__ == '__main__':
    main()
