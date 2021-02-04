import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from clustcr.chem_properties import AALPHABET
from clustcr.clustering.metrics import Metrics
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import auc, plot_roc_curve


def profile_matrix(sequences : list):
    '''
    Calculates the profile matrix for a set of sequences (i.e. all cluster members).
    NOTE: this version does not take into account the expected frequency of each amino acid at each position.
    '''

    # Amino acid alphabet
    alphabet = AALPHABET

    # Initiate profile matrix with zeros
    profile = {}
    for aa in alphabet:
        profile[aa] = [0] * len(sequences[0])
    
    # Fill in profile matrix
    for pos in range(len(sequences[0])):
        psc = pd.Series([seq[pos] for seq in sequences]).value_counts()
        for i in psc.index:
            profile[i][pos] = np.round(psc.loc[i] / len(sequences),2)
    
    # Generate output as a pd.DataFrame
    colnames = ["p" + str(p) for p in range(len(sequences[0]))]        
    profile = pd.DataFrame(profile,index=colnames).T # indices will be columns, because the df is transposed
    
    return profile


def motif_from_profile(profile):
    '''
    Generate consensus sequence motif from a profile matrix.
    Square brackets [...] indicate multiple aa possibilities at that position.
    X represents any aa.
    '''
    
    consensus = ''
    for col in profile.columns:
        if profile[col].max() > .5:
            consensus += profile[col].idxmax()
        elif sum(profile[col].nlargest(2)) >= .5:
            if profile[col].nlargest(2)[0] >= 2 * profile[col].nlargest(2)[1]:
                consensus += profile[col].idxmax()
            else:
                char = "[" + ''.join(profile[col].nlargest(2).index) + "]"
                consensus += char
        else:
            consensus += "."
    
    return consensus


def principal_component_analysis(features, labels, n_comp):
    
    # Scale features
    data = StandardScaler().fit_transform(features)
    
    # Perform PCA
    pca = PCA(n_components = n_comp)
    x_new = pca.fit_transform(data)
    
    # Plot results
    fig, ax = plt.subplots(figsize=(12,8))
    score = x_new[:,0:2]
    xs = score[:,0]
    ys = score[:,1]
    coeff = np.transpose(pca.components_[0:2, :])
    n = coeff.shape[0]
    scalex = 1.0/(xs.max() - xs.min())
    scaley = 1.0/(ys.max() - ys.min())
    plt.scatter(xs * scalex,ys * scaley, c='black')
    
    # Plot loadings
    for i in range(n):
        ax.arrow(0, 0, coeff[i,0], coeff[i,1],color = 'r',alpha = 0.5)
        ax.text(coeff[i,0]* 1.15, coeff[i,1] * 1.15, labels[i], color = 'g', ha = 'center', va = 'center', fontsize=12)
    
    # Styling
    var_expl = pca.explained_variance_ratio_
    ax.set_xlim(-1,1)
    ax.set_ylim(-1,1)
    ax.set_xlabel("PC1 ({}%)".format(np.round(var_expl[0], 3) * 100), fontsize=12)
    ax.set_ylabel("PC2 ({}%)".format(np.round(var_expl[1], 3) * 100), fontsize=12)
    ax.set_title("PCA loadings (PC 1 and PC 2)", fontsize=16)
    ax.grid()
    
    # Save figure and show output
    fig.savefig("../results/figures/cluster_features_PCA.pdf", format='pdf', dpi=1200)
    plt.show()
    
    

def generate_labels(clusters, epitope_data):
        
    metrics = Metrics(clusters, epitope_data)
    cm = metrics.calc_confmat()[0]
    
    labels = [cm[i].max() / np.sum(cm[i]) for i in cm]
    
    return labels


def data_to_ml_format(features, labels, cluster_size):
    
    X = np.array(features)
    
    y = np.asarray(labels) # Actual and permuted purities
    d = np.append(X, y, axis = 1) # Append targets to features
    d = d[d[:,1] >= cluster_size] # Filter on cluster size
    d = d[~np.isnan(d).any(axis=1)] # Remove nan values from array

    X = d[:,:-1] # Isolate features
    X = StandardScaler().fit_transform(X) # Scale features
    
    y = np.array(d[:,-1]) # Actual purities
    
    return X, y

    
def prep_data_for_ML(features, labels = None, c = .9, s = 3):
    '''
    Generate binary labels from cluster purities. These labels represent:
        - 1 ~ GOOD PURITY
        - 0 ~ BAD PURITY
    
    c : value that defines the cutoff between good and bad purity.
    s : minimum allowed cluster size.
    '''
    
    if labels is None:
        labels = generate_labels()

    for a in labels:
        # Actual clusters
        if a[1] >= c:
            labels[a[0]] = 1 # Good
        else:
            labels[a[0]] = 0 # Bad
            
    return data_to_ml_format(features, labels, s)


def stratified_cross_validation(model, X, y, n_folds = 10):
    
    # Cross-validation model
    cv = StratifiedKFold(n_splits = n_folds)
    
    # Evaluate classifier (actual purities)
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    
    fig, ax = plt.subplots(figsize=(14,11))
    for i, (train, test) in enumerate(cv.split(X, y)):
        model.fit(X[train], y[train])
        viz = plot_roc_curve(model, X[test], y[test],
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
            label=r'Mean AUC = %0.2f $\pm$ %0.2f' % (mean_auc, std_auc),
            lw=6, alpha=0.8)
    
    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='steelblue', alpha=.4,
                    label=r'$\pm$ 1 std. dev.')
    
    return fig, ax