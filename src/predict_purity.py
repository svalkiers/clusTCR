# -*- coding: utf-8 -*-
"""
author: Sebastiaan Valkiers
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor
from sklearn.decomposition import PCA
from sklearn import metrics
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import auc, plot_roc_curve


class ClusterAnalysis:

    def __init__(self, features, targets_actual, targets_permuted = None):
        
        self.feature_names = features.columns
        self.X = np.array(features)
        self.actual = np.array(targets_actual)
        self.permuted = np.array(targets_permuted)
        
        if targets_permuted is None:
            # If no permuted labels are available, permute the existing labels
            # Note that this permutation is not equal to the original one
            self.permuted = np.random.permutation(self.actual)
            
            

    def plot_purity_distribution(self, n_bins = 20, labels = ["Actual", "Permuted"]):
        '''
        Plot the distribution of cluster purities for actual and permuted clusters.
        
        n_bins : histogram bins
        labels : labels of the distributions
        '''
    
        data = np.column_stack([self.actual, self.permuted])
            
        fig, ax = plt.subplots(figsize=(12,8))
        ax.hist(data, n_bins, 
                density=True, 
                histtype='bar', 
                label=labels)
        ax.set_title("Actual versus permuted cluster purity", fontsize=16)
        ax.set_xlabel("Purity", fontsize=12)
        ax.legend(prop={'size': 14})
        plt.show()
        
        
    
    def discretize_labels(self, c = .8):
        '''
        Generate binary labels from cluster purities. These labels represent:
            - 1 ~ GOOD PURITY
            - 0 ~ BAD PURITY
        
        c : value that defines the cutoff between good and bad purity.
        '''
    
        for a, b in zip(enumerate(self.actual), enumerate(self.permuted)):
            # Actual clusters
            if a[1] >= c:
                self.actual[a[0]] = 1 # Good
            else:
                self.actual[a[0]] = 0 # Bad
            # Permuted clusters
            if b[1] >= c:
                self.permuted[b[0]] = 1 # Good
            else:
                self.permuted[b[0]] = 0 # Bad
                
        return self.actual, self.permuted
    
    
    
    def prep_data(self, s = 3):
        '''
        Prepare data for machine learning purposes.
        
        s : minimum cluster size
        '''
        
        y = np.array([self.actual, self.permuted]).T # Actual and permuted purities
        d = np.append(self.X, y, axis = 1) # Append targets to features
        d = d[d[:,1]>=s] # Filter on cluster size
        d = d[~np.isnan(d).any(axis=1)] # Remove nan values from array
    
        self.X = d[:,:-2] # Isolate features
        self.X = StandardScaler().fit_transform(self.X) # Scale features
        
        self.actual = np.array(d[:,-2]) # Actual purities
        self.permuted = np.array(d[:,-1]) # Permuted purities
        
        return self.X, self.actual, self.permuted
    
    
    
    def cluster_quality_classifier(self, n_folds = 10, feat_importances = True):
        '''
        Classification model that predicts the quality of a cluster (as defined by purity)
        based on a number of cluster features.
        
        n_folds : number of cross-validation folds (default = 10)
        feat_importances : visualize barplot with feature importances (default = True)
        '''
        
        # Define CV and RF model
        cv = StratifiedKFold(n_splits=n_folds)
        classifier = RandomForestClassifier(n_estimators=100, max_leaf_nodes=1000, criterion='entropy', 
                                            min_samples_leaf=3, min_samples_split=3, bootstrap=False, 
                                            max_depth=100, max_features='sqrt', n_jobs=-1)
        
        # Evaluate classifier (actual purities)
        tprs = []
        aucs = []
        mean_fpr = np.linspace(0, 1, 100)
        
        fig, ax = plt.subplots(figsize=(14,11))
        for i, (train, test) in enumerate(cv.split(self.X, self.actual)):
            classifier.fit(self.X[train], self.actual[train])
            viz = plot_roc_curve(classifier, self.X[test], self.actual[test],
                                 name="", alpha=0, lw=3, ax=ax, color="royalblue")
            interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
            interp_tpr[0] = 0.0
            tprs.append(interp_tpr)
            aucs.append(viz.roc_auc)
        
        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        mean_auc = auc(mean_fpr, mean_tpr)
        std_auc = np.std(aucs)
        
        # Plot ROC curve
        ax.plot(mean_fpr, mean_tpr, color='royalblue',
                label=r'Actual (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
                lw=6, alpha=0.8)
        
        std_tpr = np.std(tprs, axis=0)
        tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
        tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
        ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='steelblue', alpha=.4,
                        label=r'$\pm$ 1 std. dev.')
        
        
        # Evaluate classifier (permuted purities)
        tprs = []
        aucs = []
        mean_fpr = np.linspace(0, 1, 100)
        
        for i, (train, test) in enumerate(cv.split(self.X, self.permuted)):
            classifier.fit(self.X[train], self.permuted[train])
            viz = plot_roc_curve(classifier, self.X[test], self.permuted[test],
                                 name="", alpha=0, lw=3, ax=ax, color="royalblue")
            interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
            interp_tpr[0] = 0.0
            tprs.append(interp_tpr)
            aucs.append(viz.roc_auc)
        
        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        mean_auc = auc(mean_fpr, mean_tpr)
        std_auc = np.std(aucs)
        
        ax.plot(mean_fpr, mean_tpr, color='firebrick',
                label=r'Baseline (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
                lw=6, alpha=0.8)
        
        std_tpr = np.std(tprs, axis=0)
        tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
        tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
        ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='indianred', alpha=.2,
                        label=r'$\pm$ 1 std. dev.')
        
        # Figure aesthetics
        ax.set(xlim=[0, 1], ylim=[0, 1])
        ax.set_title("Receiver operating characteristic", fontsize=26)
        ax.set_xlabel("False Positive Rate", fontsize=20)
        ax.set_ylabel("True Positive Rate", fontsize=20)
        ax.legend(loc="lower right")
        
        plt.gcf()
        handles, labels = plt.gca().get_legend_handles_labels()
        handles = handles[10:]
        handles = handles[:1] + [handles[-3]]
        labels = labels[10:]
        labels = labels[:1] + [labels[-3]]
        by_label = dict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys(), fontsize=25)
        plt.show()
        
        if feat_importances is True:
            # Refit model on features and actual labels
            classifier.fit(self.X[train], self.actual[train])
            
            df = pd.DataFrame({"feature":self.feature_names, 
                               "importance":classifier.feature_importances_})
            df.sort_values(by="importance", ascending=False, inplace=True)
            
            fig, ax = plt.subplots(figsize=(12,8))
            sns.barplot(x=df["importance"],
                        y=df["feature"], 
                        orient="h",
                        color=sns.color_palette()[0])
            # ax.set_xticklabels(fontsize=16, labels=np.round(np.arange(0,classifier.feature_importances_.max(),0.02),2))
            ax.set_yticklabels(fontsize=16, labels=df["feature"])
            ax.set_xlabel("Importance", fontsize=20)
            ax.set_ylabel("")
            ax.set_title("Feature importance - VDJdb", fontsize=26)
    
    
    def cluster_quality_regressor(self, feat_importances = True):
        '''
        Regression model that predicts the quality of a cluster (as defined by purity)
        based on a number of cluster features.
        
        feat_importances : visualize barplot with feature importances (default = True)
        '''
        
        X_train, X_test, y_train, y_test = train_test_split(self.X, self.actual, test_size=.2)
        
        regr = RandomForestRegressor(n_estimators=100, n_jobs=-1)
        regr.fit(X_train, y_train)
        
        pred = regr.predict(X_test)
    
        print('Mean Absolute Error (MAE):', metrics.mean_absolute_error(y_test, pred))
        print('Mean Squared Error (MSE):', metrics.mean_squared_error(y_test, pred))
        print('Root Mean Squared Error (RMSE):', np.sqrt(metrics.mean_squared_error(y_test, pred)))
        mape = np.mean(np.abs((y_test - pred) / np.abs(y_test)))
        print('Mean Absolute Percentage Error (MAPE):', round(mape * 100, 2))
        print('Accuracy:', round(100*(1 - mape), 2))
    
        if feat_importances is True:
    
            fig, ax = plt.subplots(figsize=(12,8))
            df = pd.DataFrame({"feature":self.feature_names,
                               "importance":regr.feature_importances_})
            df.sort_values(by="importance", ascending=False, inplace=True)
            sns.barplot(x=df["importance"],
                        y=df["feature"], 
                        orient="h",
                        color=sns.color_palette()[0])
            # ax.set_xticklabels(fontsize=16, labels=np.round(np.arange(0,classifier.feature_importances_.max(),0.02),2))
            ax.set_yticklabels(fontsize=16, labels=df["feature"])
            ax.set_xlabel("Importance", fontsize=20)
            ax.set_ylabel("")
            ax.set_title("Feature importance - VDJdb", fontsize=26)
            plt.show()
    
    
    def perform_PCA(self, n_comp = 5):
        '''
        Perform principal component analysis using cluster features.
        
        n_comp : number of principal components (default = 5)
        '''
        
        pca = PCA(n_components = n_comp)
        x_new = pca.fit_transform(self.X)
        
        fig, ax = plt.subplots(figsize=(12,8))
        score = x_new[:,0:2]
        xs = score[:,0]
        ys = score[:,1]
        coeff = np.transpose(pca.components_[0:2, :])
        n = coeff.shape[0]
        labels = self.feature_names
        scalex = 1.0/(xs.max() - xs.min())
        scaley = 1.0/(ys.max() - ys.min())
        plt.scatter(xs * scalex,ys * scaley, c='black')
        for i in range(n):
            ax.arrow(0, 0, coeff[i,0], coeff[i,1],color = 'r',alpha = 0.5)
            ax.text(coeff[i,0]* 1.15, coeff[i,1] * 1.15, labels[i], color = 'g', ha = 'center', va = 'center', fontsize=12)
        ax.set_xlim(-1,1)
        ax.set_ylim(-1,1)
        ax.set_xlabel("PC{}".format(1), fontsize=12)
        ax.set_ylabel("PC{}".format(2), fontsize=12)
        ax.set_title("PCA loadings (PC 1 and PC 2)", fontsize=16)
        ax.grid()
        plt.show()