from base import TestBase
from clustcr import datasets, read_cdr3, metarepertoire


class InputTest(TestBase):

    def setUp(self):
        pass

    def test_testing_data(self):
        datasets.test_cdr3()
        datasets.test_epitopes()

    def test_vdj(self):
        datasets.vdjdb_cdr3()
        datasets.vdjdb_epitopes()

    def test_read_immune(self):
        read_cdr3('input/immuneaccess/HIP05763.tsv', data_format='immuneaccess')

    def test_metarepertoire(self):
        metarepertoire('input/immuneaccess', data_format='immuneaccess', out_format='CDR3', n_sequences=10 ** 3)

