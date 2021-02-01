import numpy as np 
import pickle 
import matplotlib.pyplot as plt 
 
from clustcr.analysis.tools import principal_component_analysis, stratified_cross_validation 
from clustcr.analysis.feature_generator import Features 
from sklearn.preprocessing import StandardScaler 
from sklearn.ensemble import RandomForestClassifier 
 
 
class cluster_analysis: 
     
    def __init__(self, features): 
        self.features = features 
         
     
    def _pca(self, number_of_components = 5): 
        ''' 
        Perform principal component analysis using cluster features. 
        ''' 
         
        labels = self.features.columns 
        X = np.array(self.features) 
        X = X[~np.isnan(X).any(axis=1)] 
         
        principal_component_analysis(X, labels, number_of_components) 
         
     
    def _predict_quality(features : np.array, model = None): 
        """ 
        Predict clustering quality from a set of clustering features. 
        A pre-trained, default model is provides, but customization is possible. 
         
        Parameters 
        ---------- 
        features: np.array 
            Numpy array containing clustering features. These can be calculated 
            using the feature_generator module. 
        model: .pkl file, optional 
            Prediction model. Default model is provided. Custom models with a 
            .pkl extension can also be used. 
             
        Returns 
        ------- 
        predictions: np.array 
            Returns a numpy array indicating high (1) or low (0) quality 
            of clusters. 
        """ 
         
        if model is None: 
            with open('./cq_classifier.pkl', 'rb') as f: 
                model = pickle.load(f) 
             
        X = StandardScaler().fit_transform(features) 
        predictions = model.predict(X) 
         
        return predictions 
     
 
class train_model: 
     
    def __init__(self, clusters, epitope_data, features=None): 
         
        self.clusters = clusters 
        self.epitopes = epitope_data 
         
        if features is not None: 
            self.features = features 
        else: 
            f = Features(self.clusters) 
            self.features = f.compute_features() 
         
         
    def fit(self): 
     
        data = self.prep_data_for_ML() 
         
        X = data[0] 
        y = data[1] 
         
        # Define CV and RF model 
        classifier = RandomForestClassifier(n_estimators=100, max_leaf_nodes=1000, criterion='entropy',  
                                            min_samples_leaf=3, min_samples_split=3, bootstrap=False,  
                                            max_depth=100, max_features='sqrt', n_jobs=-1) 
         
        return classifier.fit(X, y) 
     
     
    def evaluate(self): 
         
        data = self.prep_data_for_ML() 
         
        X = data[0] 
        y = data[1] 
         
        classifier = RandomForestClassifier(n_estimators=100, max_leaf_nodes=1000, criterion='entropy',  
                                            min_samples_leaf=3, min_samples_split=3, bootstrap=False,  
                                            max_depth=100, max_features='sqrt', n_jobs=-1) 
         
        fig, ax = stratified_cross_validation(classifier, X, y) 
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