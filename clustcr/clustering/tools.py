import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import auc, plot_roc_curve

from clustcr.chem_properties import AALPHABET


def create_edgelist(cdr3):
    '''
    Create tab-separated edgelist of edges with HD = 1, from a set of sequences.    
    '''
    # Set makes sure there are no dupes
    cdr3 = set(cdr3)
    
    # Hashing
    cdr3hash = dict()
    for cdr in cdr3:
        for hash in (cdr[::2], cdr[1::2]):
            if hash not in cdr3hash:
                cdr3hash[hash] = set()
            cdr3hash[hash].add(cdr)
            
    # Generate network
    edgelist = set()
    for hash in cdr3hash:
        if len(cdr3hash[hash]) >= 1:
            for cdr1 in cdr3hash[hash]:
                for cdr2 in cdr3hash[hash]:
                    if cdr1 != cdr2:
                        if cdr1 <= cdr2:
                            if sum(ch1 != ch2 for ch1, ch2 in zip(cdr1, cdr2)) <= 1:
                                edgelist.add(cdr1 + "\t" + cdr2)

    return edgelist


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
            consensus += "X"
    
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
    for i in range(n):
        ax.arrow(0, 0, coeff[i,0], coeff[i,1],color = 'r',alpha = 0.5)
        ax.text(coeff[i,0]* 1.15, coeff[i,1] * 1.15, labels[i], color = 'g', ha = 'center', va = 'center', fontsize=12)
    var_expl = pca.explained_variance_ratio_
    ax.set_xlim(-1,1)
    ax.set_ylim(-1,1)
    ax.set_xlabel("PC1 ({}%)".format(np.round(var_expl[0], 3) * 100), fontsize=12)
    ax.set_ylabel("PC2 ({}%)".format(np.round(var_expl[1], 3) * 100), fontsize=12)
    ax.set_title("PCA loadings (PC 1 and PC 2)", fontsize=16)
    ax.grid()
    fig.savefig("../results/figures/cluster_features_PCA.pdf", format='pdf', dpi=1200)
    plt.show()
    
    
def data_to_ml_format(features, actual_labels, permuted_labels, cluster_size):
    
    X = np.array(features)
    
    y = np.array([actual_labels, permuted_labels]).T # Actual and permuted purities
    d = np.append(X, y, axis = 1) # Append targets to features
    d = d[d[:,1] >= cluster_size] # Filter on cluster size
    d = d[~np.isnan(d).any(axis=1)] # Remove nan values from array

    X = d[:,:-2] # Isolate features
    X = StandardScaler().fit_transform(X) # Scale features
    
    y_a = np.array(d[:,-2]) # Actual purities
    y_p = np.array(d[:,-1]) # Permuted purities
    
    return X, y_a, y_p


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
