import pandas as pd


class MultiRepertoireFeatureMatrix:

    def __init__(self):
        self.feature_dict = {}

    def get_matrix(self):
        df = pd.DataFrame.from_dict(self.feature_dict, orient='index')
        return df.fillna(0).reindex(sorted(df.columns), axis=1).astype(int)

    def add(self, preclusters, clusters):
        name_dict = self._make_name_dict(preclusters)
        for index, row in clusters.iterrows():
            cdr3, cluster = row['CDR3'], row['cluster']
            for name in name_dict.get(cdr3, []):
                self._add_to_feature_dict(name, cluster)

    def _make_name_dict(self, preclusters):
        name_dict = {}
        for index, row in preclusters.clusters_df.iterrows():
            cdr3, name = row['CDR3'], row['name']
            if cdr3 not in name_dict:
                name_dict[cdr3] = set()
            name_dict[cdr3].add(name)
        return name_dict

    def _add_to_feature_dict(self, name, cluster):
        if name not in self.feature_dict:
            self.feature_dict[name] = {}
        entry = self.feature_dict[name]
        if cluster not in entry:
            entry[cluster] = 1
        else:
            entry[cluster] += 1

