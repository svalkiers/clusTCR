import pandas as pd
import numpy as np


class Metrics:
    """
    Metrics for evaluating clustering quality
    You can only calculate these cluster metrics if the CDR3 sequences are labelled (i.e. epitope specificity is known)!
    """

    def __init__(self, nodelist, epidata, name=None):
        self.nodelist = nodelist  # pd.DataFrame with columns ["CDR3", "cluster"]
        self.epidata = epidata  # 'ground truth', pd.DataFrame of CDR3 sequences with corresponding epitope specificity (columns=["CDR3", "Epitope"])

        # Ensure all values correspond to CDR3s in nodelist and no duplicates remain
        self.gt = self.epidata[self.epidata["CDR3"].isin(self.nodelist["CDR3"])]
        self.gt = self.gt.drop_duplicates()

        # Construct joint pd.DataFrame that stores information about cluster and epitope association of CDR3s
        self.gt = pd.merge(left=self.epidata, right=self.nodelist, on="CDR3")

        # Make a copy and permute cluster assignment column, this provides a baseline for comparison
        self.gt_baseline = self.gt.copy()
        self.gt_baseline["cluster"] = np.random.permutation(self.gt_baseline["cluster"])

        self.name = name

    def calc_confmat(self):
        """
        Construct confusion matrices for true and baseline.
        """
        self.gt["count"] = 1
        self.gt_baseline["count"] = 1
        conf_mat_t = pd.pivot_table(self.gt, values='count',
                                    index=self.gt["Epitope"],
                                    columns=self.gt["cluster"],
                                    aggfunc=np.sum,
                                    fill_value=0)
        conf_mat_b = pd.pivot_table(self.gt_baseline,
                                    values='count',
                                    index=self.gt_baseline["Epitope"],
                                    columns=self.gt_baseline["cluster"],
                                    aggfunc=np.sum,
                                    fill_value=0)

        return conf_mat_t, conf_mat_b

    def retention(self):
        '''
        Cluster retention is the fraction of sequences that has been assigned to any cluster.
        '''
        return len(self.nodelist.CDR3.unique()) / len(self.epidata.CDR3.unique())

    def purity(self, conf_mat=None):
        '''
        Method that estimates the precision of the solution.
        We assigned each cluster to the most common epitope.
        All other epitopes in the same cluster are considered false positives.
        '''

        if conf_mat is None:
            conf_mat = self.calc_confmat()

        hits_t = np.sum(conf_mat[0].apply(np.max, axis=0))
        hits_b = np.sum(conf_mat[1].apply(np.max, axis=0))

        return (hits_t / np.sum(conf_mat[0].values, axis=None), hits_b / np.sum(conf_mat[1].values, axis=None))

    def purity_90(self, conf_mat=None):
        """
        Method that determines the fraction of clusters that have
        a purity greater than .90
        """
        if conf_mat is None:
            conf_mat = self.calc_confmat()

        c_t = 0
        for col in conf_mat[0].columns:
            p_t = conf_mat[0][col].max() / conf_mat[0][col].sum()
            if p_t >= .9:
                c_t += 1
        c_b = 0
        for col in conf_mat[1].columns:
            p_b = conf_mat[1][col].max() / conf_mat[1][col].sum()
            if p_b >= .9:
                c_b += 1

        return (c_t / len(conf_mat[0].columns), c_b / len(conf_mat[1].columns))

    def consistency(self, conf_mat=None):
        """
        Method that pretends that we solved a supervised problem where each cluster corresponds to a single epitope.
        Returns the accuracy of the best solution.
        """

        if conf_mat is None:
            conf_mat = self.calc_confmat()

        # Define recursive function that finds the best fit for the diagonal
        def rec_max(mat):
            high = mat.max().max()
            col = mat.max().idxmax()
            row = mat[col].idxmax()

            if (len(mat.index) > 1 and len(mat.columns) > 1):
                high = high + rec_max(mat.drop(row, axis=0).drop(col, axis=1))

            return high

        return (rec_max(conf_mat[0]) / len(self.gt), rec_max(conf_mat[1]) / len(self.gt_baseline))

    def summary(self, conf_mat=None):
        """
        Calculates all available clustering metrics and outputs them as
        a pd.DataFrame.

        Parameters
        ----------
        conf_mat : np.array, optional
            Confusion matrix. This is automatically calculated if not provided.

        Returns
        -------
        summ : pd.DataFrame
            Summary of clustering metrics.

        """

        summ = pd.DataFrame({'actual': [self.retention(),
                                        self.purity(conf_mat)[0],
                                        self.purity_90(conf_mat)[0],
                                        self.consistency(conf_mat)[0]],
                             'baseline': [self.retention(),
                                          self.purity(conf_mat)[1],
                                          self.purity_90(conf_mat)[1],
                                          self.consistency(conf_mat)[1]], })
        summ['metrics'] = ['retention', 'purity', 'purity_90', 'consistency']

        if self.name is not None:
            summ['method'] = [self.name] * len(summ)

        return summ

