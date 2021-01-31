from clustcr import datasets

cdr3 = datasets.metarepertoire_cdr3('/home/max/PycharmProjects/clusTCR/data/immuneACCESS', 100 * (10 ** 6))
print(len(cdr3))
cdr3.to_csv()

