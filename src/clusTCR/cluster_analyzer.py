import numpy as np
import matplotlib.pyplot as plt
import pickle

from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier

from src.clusTCR.clustering.clusTCR import Features, Metrics
from toolkit.tools import principal_component_analysis, data_to_ml_format, stratified_cross_validation

class FeatureGenerator:
    
    def __init__(self, cluster_output):
        self.cluster_output = cluster_output
        
        
        
    def get_features(self):
        
        features = Features(self.cluster_output)
        aavar = features.calc_variation()
        pchem = features.calc_physchem()
        pgen = features.calc_pgen()
        
        output = features.combine(aavar, pchem, pgen)
        
        return output
    
    
    
    def perform_PCA(self, cluster_features, number_of_components = 5):
        '''
        Perform principal component analysis using cluster features.
        '''
        
        labels = cluster_features.columns
        X = np.array(cluster_features)
        X = X[~np.isnan(X).any(axis=1)]
        
        principal_component_analysis(X, labels, number_of_components)
        


class QualityPredict:
    
    def __init__(self, cluster_features):
        self.features = cluster_features
        
    
    
    def predict_clustering_quality(self, model = None):
        
        if model is None:
            with open('toolkit/cq_classifier.pkl', 'rb') as f:
                model = pickle.load(f)
            
        X = StandardScaler().fit_transform(self.features)
        predictions = model.predict(X)
        
        return predictions



class ModelTraining:
    
    def __init__(self, cluster_features, cluster_output, epitope_data):
        self.features = cluster_features
        self.clusters = cluster_output
        self.epidata = epitope_data
        
        
    
    def generate_labels(self):
        
        metrics = Metrics(self.clusters, self.epidata)
        cm_actual, cm_permuted = metrics.calc_confmat()
        
        actual = [cm_actual[i].max() / np.sum(cm_actual[i]) for i in cm_actual]
        permuted = [cm_permuted[i].max() / np.sum(cm_permuted[i]) for i in cm_permuted]
        
        return actual, permuted
    
    
    
    def plot_purity_distribution(self, actual = None, permuted = None):
        '''
        Plot the distribution of cluster purities for actual and permuted clusters.
        '''
        
        if actual is None:
            actual, permuted = self.generate_labels()
            
        data = np.column_stack([actual, permuted])
        labels = ["Actual", "Permuted"]
        
        fig, ax = plt.subplots(figsize=(12,8))
        ax.hist(data, bins = 20, 
                density = True, 
                histtype = 'bar', 
                label = labels,
                color = ['steelblue', 'lightcoral'])
        ax.set_title("Actual versus permuted cluster purity", fontsize=16)
        ax.set_xlabel("Purity", fontsize=12)
        ax.legend(prop={'size': 14})
        plt.savefig("../results/figures/purity_distribution.pdf", format='png', dpi=1200)
        plt.show()
        
        
        
    def prep_data_for_ML(self, actual = None, permuted = None, c = .9, s = 3):
        '''
        Generate binary labels from cluster purities. These labels represent:
            - 1 ~ GOOD PURITY
            - 0 ~ BAD PURITY
        
        c : value that defines the cutoff between good and bad purity.
        s : minimum allowed cluster size.
        '''
        
        if actual is None:
            actual, permuted = self.generate_labels()
    
        for a, b in zip(enumerate(actual), enumerate(permuted)):
            # Actual clusters
            if a[1] >= c:
                actual[a[0]] = 1 # Good
            else:
                actual[a[0]] = 0 # Bad
            # Permuted clusters
            if b[1] >= c:
                permuted[b[0]] = 1 # Good
            else:
                permuted[b[0]] = 0 # Bad

        return data_to_ml_format(self.features, actual, permuted, s)
    
    
    
    def train_model(self):
        
        data = self.prep_data_for_ML()
        
        X = data[0]
        y = data[1]
        
        # Define CV and RF model
        classifier = RandomForestClassifier(n_estimators=100, max_leaf_nodes=1000, criterion='entropy', 
                                            min_samples_leaf=3, min_samples_split=3, bootstrap=False, 
                                            max_depth=100, max_features='sqrt', n_jobs=-1)
        
        return classifier.fit(X, y)
    
    
    
    def evaluate_model(self):
        
        data = self.prep_data_for_ML()
        
        X = data[0]
        actual = data[1]
        # permuted = data[2]
        
        classifier = RandomForestClassifier(n_estimators=100, max_leaf_nodes=1000, criterion='entropy', 
                                            min_samples_leaf=3, min_samples_split=3, bootstrap=False, 
                                            max_depth=100, max_features='sqrt', n_jobs=-1)
        
        fig, ax = stratified_cross_validation(classifier, X, actual)
        ax.plot([0, 1], [0, 1], color = 'red', alpha = .5,  lw = 4, linestyle = '--', label = 'random')
        
        # Figure aesthetics
        ax.set(xlim=[0, 1], ylim=[0, 1])
        ax.set_title("Receiver operating characteristic", fontsize=26)
        ax.set_xlabel("False Positive Rate", fontsize=20)
        ax.set_ylabel("True Positive Rate", fontsize=20)
        ax.legend(loc="lower right")
        
        plt.gcf()
        handles, labels = plt.gca().get_legend_handles_labels()
        handles = handles[10:]
        labels = labels[10:]
        by_label = dict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys(), fontsize=25)
        fig.savefig("../results/figures/cluster_quality_ROC.pdf", format='pdf', dpi=1200)
        plt.show()