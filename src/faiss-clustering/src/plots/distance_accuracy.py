import matplotlib.pyplot as plt
import numpy as np
from src.distance import createNetwork, make_pairs_with_distance, distance1_accuracy, make_pairs_with_distance_after_sorting


def plot_accuracy(data, dataset_name):
    plt.figure()
    total_pairs_of_distance1 = len(createNetwork(data))

    x = []
    y_sorted_faiss = []
    y_faiss = []

    for items_per_cluster in np.arange(10, 110, 30):
        x.append(items_per_cluster)

        pairs_sorted_faiss, contents = make_pairs_with_distance_after_sorting(data, int(items_per_cluster))
        accuracy = distance1_accuracy(pairs_sorted_faiss, data, total_pairs_of_distance1)
        y_sorted_faiss.append(accuracy)

        pairs_faiss, contents = make_pairs_with_distance(data, int(items_per_cluster))
        accuracy = distance1_accuracy(pairs_faiss, data, total_pairs_of_distance1)
        y_faiss.append(accuracy)

        # Progress
        print(items_per_cluster)

    plt.plot(x, y_faiss, label='Faiss method')
    plt.plot(x, y_sorted_faiss, label='Faiss method (after sorting on length)')

    plt.xlabel('Average amount of items per cluster')
    plt.ylabel('Percentage of pairs found')
    plt.ylim(0, 1)
    plt.legend(loc='best')
    plt.title(f'Finding pairs with distance 1 (in {dataset_name})')
    plt.show()