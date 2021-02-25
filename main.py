from clustcr import datasets, Clustering
from clustcr.clustering.clustering import ClusteringResult
from clustcr.modules.gliph2 import gliph2
from clustcr.modules.ismart import ismart
import matplotlib.pyplot as plt
import time


def plot(metric_lambda, metric_name):
    plt.figure()
    x = []

    ismart_y, gliph2_y, ismart_clustcr_y, gliph2_clustcr_y, mcl_y = [], [], [], [], []
    for q in [0, 1, 2]:
        epitopes = datasets.vdjdb_epitopes_small(q=q)
        data_gliph2 = datasets.vdjdb_gliph2_small(q=q)
        data = datasets.vdjdb_cdr3_small(q=q)
        x.append(len(data))

        ismart_clustcr_y.append(metric_lambda(Clustering(second_step='ismart', n_cpus='all').fit(data).metrics(epitopes)))
        gliph2_clustcr_y.append(metric_lambda(Clustering(second_step='gliph2', n_cpus='all').fit(data_gliph2).metrics(epitopes)))
        mcl_y.append(metric_lambda(Clustering(second_step='mcl').fit(data).metrics(epitopes)))
        ismart_y.append(metric_lambda(ClusteringResult(ismart.iSMART(data)[0]).metrics(epitopes)))
        gliph2_y.append(metric_lambda(ClusteringResult(gliph2.GLIPH2(data_gliph2)[0]).metrics(epitopes)))

    plt.plot(x, ismart_y, label='iSMART', marker='o', linewidth=4)
    plt.plot(x, gliph2_y, label='GLIPH2', marker='o', linewidth=4)
    plt.plot(x, ismart_clustcr_y, label='Two Step with iSMART', marker='o', linewidth=2)
    plt.plot(x, gliph2_clustcr_y, label='Two Step with GLIPH2', marker='o', linewidth=2)
    plt.plot(x, mcl_y, label='Two Step with MCL', marker='o', linewidth=2)
    plt.xlabel('Number of input sequences')
    plt.ylabel(metric_name)
    plt.title('Clustering VDJdb - method comparison')
    plt.legend(loc='best')
    plt.show()


def timer(func):
    start = time.time()
    func()
    return time.time() - start


def plot_time():
    plt.figure()
    x = []

    ismart_y, gliph2_y, ismart_clustcr_y, ismart_par_clustcr_y, gliph2_clustcr_y, gliph2_par_clustcr_y, mcl_y = [], [], [], [], [], [], []
    for q in [500, 10000, 100000]:
        # data_gliph2 = datasets.vdjdb_gliph2_small(q=q)
        # data = datasets.vdjdb_cdr3_small(q=q)
        data_gliph2 = datasets.metarepertoire('/home/max/Downloads/emerson-2017-natgen', 'immuneaccess', out_format='GLIPH2', n_sequences=q)
        data = datasets.metarepertoire('/home/max/Downloads/emerson-2017-natgen', 'immuneaccess', n_sequences=q)
        x.append(len(data))

        ismart_clustcr_y.append(
            timer(lambda: Clustering(second_step='ismart').fit(data)))
        gliph2_clustcr_y.append(
            timer(lambda: Clustering(second_step='gliph2').fit(data_gliph2)))
        ismart_par_clustcr_y.append(
            timer(lambda: Clustering(second_step='ismart', n_cpus='all').fit(data)))
        gliph2_par_clustcr_y.append(
            timer(lambda: Clustering(second_step='gliph2', n_cpus='all').fit(data_gliph2)))
        mcl_y.append(timer(lambda: Clustering(second_step='mcl').fit(data)))
        ismart_y.append(timer(lambda: ClusteringResult(ismart.iSMART(data)[0])))
        gliph2_y.append(timer(lambda: ClusteringResult(gliph2.GLIPH2(data_gliph2)[0])))

    plt.plot(x, ismart_y, label='iSMART', marker='o')
    plt.plot(x, gliph2_y, label='GLIPH2', marker='o')
    plt.plot(x, ismart_clustcr_y, label='Two Step iSMART', marker='o')
    plt.plot(x, gliph2_clustcr_y, label='Two Step GLIPH2', marker='o')
    plt.plot(x, ismart_par_clustcr_y, label='Two Step iSMART in parallel', marker='o')
    plt.plot(x, gliph2_par_clustcr_y, label='Two Step GLIPH2 in parallel', marker='o')
    plt.plot(x, mcl_y, label='Two Step MCL', marker='o')
    plt.xlabel('Number of input sequences')
    plt.ylabel('Time (in seconds)')
    plt.title('Clustering VDJdb - method comparison')
    plt.legend(loc='best')
    plt.show()


plot_time()
plot(lambda x: x.purity()[0], 'Purity')
plot(lambda x: x.consistency()[0], 'Consistency')
plot(lambda x: x.retention(), 'Retention')

