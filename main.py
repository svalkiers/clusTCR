from clustcr import Clustering, datasets
import time
import matplotlib.pyplot as plt

DATA = []
for q in [0, 1]:
    DATA.append((datasets.vdjdb_cdr3_small(q), datasets.vdjdb_epitopes_small(q)))


def plot(metric_lambda, name):
    plt.figure()
    x = []
    y = {}
    y_standard = []
    for data, epi in DATA:
        x.append(len(data))
        starting_time = time.time()
        clust = Clustering(n_cpus='all', blosum_cutoff=None).fit(data)
        y_standard.append(time.time() - starting_time)
        for bl in [0.9, 0.93, 0.96]:
            for hd in [1, 2, 3]:
                starting_time = time.time()
                clustering = Clustering(n_cpus='all', blosum_cutoff=bl, hd_cutoff=hd).fit(data)
                key = (bl, hd)
                if key not in y:
                    y[key] = []
                y[key].append(time.time() - starting_time)

    for key in y:
        linestyle = ['dotted', 'dashed', 'dashdot'][int(key[1]) - 1]
        plt.plot(x, y[key], label=f'HD {key[1]}, BL {key[0]}', linestyle=linestyle)
    plt.plot(x, y_standard, label='HD 1, BL 0')
    plt.xlabel('Number of input sequences')
    plt.ylabel(name)
    plt.legend(loc='best')
    plt.title(f'Hyperparameter tuning')
    plt.show()


def plotHDtime_parallel():
    plt.figure()
    x = []
    y = {}
    for data, epi in DATA:
        x.append(len(data))
        for hd_multiprocessing in [False, True]:
            for hd in [1, 2, 3]:
                print(hd)
                clustering = Clustering(n_cpus='all', blosum_cutoff=None, hd_cutoff=hd)
                clustering.hd_multiprocessing = hd_multiprocessing
                starting_time = time.time()
                clustering.fit(data)
                duration = time.time() - starting_time
                key = (hd, hd_multiprocessing)
                if key not in y:
                    y[key] = []
                y[key].append(duration)

    for key in y:
        linestyle = ['dotted', 'dashed', 'dashdot'][int(key[0]) - 1]
        plt.plot(x, y[key], label=f'HD {key[0]} {"parallel" if key[1] else ""}', linestyle=linestyle)
    plt.xlabel('Number of input sequences')
    plt.ylabel('Time')
    plt.legend(loc='best')
    plt.title(f'Hyperparameter tuning')
    plt.show()


if __name__ == '__main__':
    # plot(lambda x, y: x.purity()[0], 'Purity')
    # plot(lambda x, y: x.consistency()[0], 'Consistency')
    # plot(lambda x, y: x.retention(), 'Retention')
    # plot(lambda x, y: time.time() - y, 'Time')
    plotHDtime_parallel()
