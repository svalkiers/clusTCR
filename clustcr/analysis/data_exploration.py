import numpy as np 
import pickle 
import matplotlib.pyplot as plt 
 
from clustcr.analysis.tools import principal_component_analysis, stratified_cross_validation 
from clustcr.analysis.features import FeatureGenerator
from clustcr.clustering.metrics import Metrics
from sklearn.preprocessing import StandardScaler 
from sklearn.ensemble import RandomForestClassifier 
from os.path import join, dirname, abspath

DIR = dirname(abspath(__file__)) 
 
class ClusterAnalysis:
    """
    Analysis module that can be used for the exploration of clustering results.
    Analyses are based on cluster features, calculated using the FeatureGenerator.
    ClusterAnalyses provides methods for performing PCA, and predicting
    clustering quality.
    """
     
    def __init__(self, features): 
        self.features = features 
         
     
    def pca(self, number_of_components=5, name=None): 
        """ 
        Perform principal component analysis using cluster features. 
        
        Parameters
        ----------
        number_of_components: int
            Map features to n number of principal components.
            Default is 5.
        """ 
         
        labels = self.features.columns 
        X = np.array(self.features) 
        X = X[~np.isnan(X).any(axis=1)] 
         
        principal_component_analysis(X, labels, number_of_components, name)
         
     
    def predict_quality(self, model=None): 
        """ 
        Predict clustering quality from a set of clustering features. 
        A pre-trained, default model is provided, but customization is possible.
        To use a custom model, provide the path to the model with a .pkl
        extension. The model can be created using the TrainModel class.
         
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
        
        assert 'pgen_avg' in self.features.columns, 'Error: classifier requires values for pgen_avg and pgen_var.'
        
        if model is None: 
            with open(join(DIR,'cq_classifier.pkl'),'rb') as f: 
                model = pickle.load(f)
        else:
            with open(model, 'rb') as f:
                model = pickle.load(f)
             
        X = StandardScaler().fit_transform(self.features) 
        predictions = model.predict(X) 
         
        return predictions 
     
 
class TrainModel:
    """
    Train your own classification model for predicting clustering quality.
    This requires information about cluster, epitopes and the features
    corresponding to the clusters. If the latter is not provided, clusTCR
    will calculate the features of the clustering data provided.
    """
     
    def __init__(self, clusters, epitope_data, features=None): 
         
        self.clusters = clusters 
        self.epitopes = epitope_data 
         
        if features is not None: 
            self.features = features 
        else: 
            f = FeatureGenerator(self.clusters) 
            self.features = f.compute_features()
    
    
    def _generate_labels(self):
        
        metrics = Metrics(self.clusters, self.epitopes)
        cm = metrics.calc_confmat()[0]
        
        labels = [cm[i].max() / np.sum(cm[i]) for i in cm]
        
        return labels
            
            
    def _data_to_ml_format(self, labels, cluster_size):
    
        d = self.features
        y = np.asarray(labels)
        
        d['y'] = labels
        d = np.asarray(d)
        # y = np.asarray(labels) # Actual and permuted purities
        # d = np.append(X, y, axis = 1) # Append targets to features
        d = d[d[:,1] >= cluster_size] # Filter on cluster size
        d = d[~np.isnan(d).any(axis=1)] # Remove nan values from array
    
        X = d[:,:-1] # Isolate features
        X = StandardScaler().fit_transform(X) # Scale features
        
        y = np.array(d[:,-1]) # Actual purities
        
        return X, y


    def _prep_data_for_ML(self, labels = None, c = .9, s = 3):
        '''
        Generate binary labels from cluster purities. These labels represent:
            - 1 ~ GOOD PURITY
            - 0 ~ BAD PURITY
        
        c : value that defines the cutoff between good and bad purity.
        s : minimum allowed cluster size.
        '''
        
        if labels is None:
            labels = self._generate_labels()
    
        for i,j in enumerate(labels):
            # Actual clusters
            if j >= c:
                labels[i] = 1 # Good
            else:
                labels[i] = 0 # Bad
                
        return self._data_to_ml_format(labels, s)

         
    def fit_data(self):
        """
        Fit your model to the data.
        """
     
        data = self._prep_data_for_ML() 
         
        X = data[0] 
        y = data[1] 
         
        # Define CV and RF model 
        classifier = RandomForestClassifier(n_estimators=100, max_leaf_nodes=1000, criterion='entropy',  
                                            min_samples_leaf=3, min_samples_split=3, bootstrap=False,  
                                            max_depth=100, max_features='sqrt', n_jobs=-1) 
         
        return classifier.fit(X, y), X, y
     
     
    def evaluate(self, location=None):
        """
        Evaluate your model through stratified cross-validation procedure.
        Model performance is expressed as the area under the receiver
        operating characteristic (AUROC).
        """
         
        data = self._prep_data_for_ML() 
         
        X = data[0] 
        y = data[1] 
         
        classifier = RandomForestClassifier(n_estimators=100, max_leaf_nodes=1000, criterion='entropy',  
                                            min_samples_leaf=3, min_samples_split=3, bootstrap=False,  
                                            max_depth=100, max_features='sqrt', n_jobs=-1) 
         
        fig, ax = stratified_cross_validation(classifier, X, y) 
        ax.plot([0, 1], [0, 1], color = 'red', alpha = .5,  lw = 4, linestyle = '--', label = 'random') 
         
        # Plot styling 
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
        plt.show()
        
        if location is not None:
            fig.savefig(location + 'cq_predict_roc.eps', format='eps')

    
    def save(self, model, filename):
        """
        Dump custom model into .pkl file.
        
        Parameters
        ----------
        model:  
            Classification model trained with the ._fit() method.
        filename: str
            Path to location where model needs to be dumped.
        """
        pickle.dump(model, open(filename, 'wb'))
