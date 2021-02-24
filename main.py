from clustcr import datasets, Clustering
import matplotlib.pyplot as plt

ismart = Clustering(second_step='ismart')
gliph2 = Clustering(second_step='gliph2')
mcl = Clustering(second_step='mcl')


def plot(metric_lambda, metric_name):
    plt.figure()
    x = []
    ismart_y = []
    gliph2_y = []
    mcl_y = []
    for q in [0, 1, 2]:
        epitopes = datasets.vdjdb_epitopes_small(q=q)
        data_gliph2 = datasets.vdjdb_gliph2_small(q=q)
        data = datasets.vdjdb_cdr3_small(q=q)
        x.append(len(data))

        ismart_metrics = ismart.fit(data).metrics(epitopes)
        ismart_y.append(metric_lambda(ismart_metrics))

        gliph2_metrics = gliph2.fit(data_gliph2).metrics(epitopes)
        gliph2_y.append(metric_lambda(gliph2_metrics))

        mcl_metrics = mcl.fit(data).metrics(epitopes)
        mcl_y.append(metric_lambda(mcl_metrics))

    plt.plot(x, ismart_y, label='iSMART', marker='o')
    plt.plot(x, gliph2_y, label='GLIPH2', marker='o')
    plt.plot(x, mcl_y, label='MCL', marker='o')
    plt.xlabel('Number of input sequences')
    plt.ylabel(metric_name)
    plt.title('clusTCR second step - method comparison')
    plt.legend(loc='best')
    plt.show()


plot(lambda x: x.purity()[0], 'Purity')
plot(lambda x: x.consistency()[0], 'Consistency')
plot(lambda x: x.retention(), 'Retention')




