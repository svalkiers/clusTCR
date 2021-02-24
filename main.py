from clustcr import datasets, Clustering
from clustcr.clustering.clustering import ClusteringResult
from clustcr.modules.gliph2 import gliph2
from clustcr.modules.ismart import ismart
import matplotlib.pyplot as plt

ismart_clustcr = Clustering(second_step='ismart')
gliph2_clustcr = Clustering(second_step='gliph2')
mcl = Clustering(second_step='mcl')


def plot(metric_lambda, metric_name):
    plt.figure()
    x = []

    ismart_y, gliph2_y, ismart_clustcr_y, gliph2_clustcr_y, mcl_y = [], [], [], [], []
    for q in [0, 1, 2]:
        epitopes = datasets.vdjdb_epitopes_small(q=q)
        data_gliph2 = datasets.vdjdb_gliph2_small(q=q)
        data = datasets.vdjdb_cdr3_small(q=q)
        x.append(len(data))

        ismart_clustcr_y.append(metric_lambda(ismart_clustcr.fit(data).metrics(epitopes)))
        gliph2_clustcr_y.append(metric_lambda(gliph2_clustcr.fit(data_gliph2).metrics(epitopes)))
        mcl_y.append(metric_lambda(mcl.fit(data).metrics(epitopes)))
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


plot(lambda x: x.purity()[0], 'Purity')
plot(lambda x: x.consistency()[0], 'Consistency')
plot(lambda x: x.retention(), 'Retention')




