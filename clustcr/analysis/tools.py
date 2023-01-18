import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from collections import Counter

from clustcr.chem_properties import AALPHABET
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import auc, roc_curve


def profile_matrix(sequences : list):
    '''
    Calculates the profile matrix for a set of sequences (i.e. all cluster members).
    NOTE: this version does not take into account the expected frequency of each amino acid at each position.
    '''

    # Make sure to proceed only if all sequences in the cluster have equal length
    seq_len = len(sequences[0])
    if not all(len(seq) == seq_len for seq in sequences):

        # On the rare occasion that a cluster contains sequences of inequal length.
        # Typically, there is/are only one (or very few) sequence(s) that differ from the avg. CDR3 length in the cluster.
        # Therefore, we use the length of the highest proportion of sequences as the standard, and delete all others.
        seq_len = Counter([len(s) for s in sequences]).most_common()[0][0]
        sequences = [s for s in sequences if len(s) == seq_len]
    
    # Initiate profile matrix with zeros
    pm = np.zeros(shape=(len(AALPHABET), seq_len))

    # initiate AA dict:
    AAs = {aa: i for i, aa in enumerate(AALPHABET)}

    # Fill in profile matrix with counts
    for s in sequences:
        for i, aa in enumerate(s):
            pm[AAs[aa], i] += 1

    # normalize profile matrix to percentages
    pm = pm / len(sequences)

    return pm


def motif_from_profile(profile, method, cutoff=.7):
    '''
    Generate consensus sequence motif from a profile matrix.
    Square brackets [...] indicate multiple aa possibilities at that position.
    X represents any aa.
    '''
    AA_map = {i:aa for i,aa in enumerate(AALPHABET)}

    consensus = ''
    
    if method.lower() == 'standard':
        top_idxs = np.argpartition(profile, -2, axis=0)[-2:].T
        top_values = np.partition(profile, -2, axis=0)[-2:].T
        for (second_max_idx, max_idx), (second_max_value, max_value) in zip(top_idxs, top_values):
            if max_value >= cutoff:
                consensus += AA_map[max_idx]
            elif max_value + second_max_value >= cutoff:
                if max_value >= 2*second_max_value:
                    consensus += AA_map[max_idx].lower()
                else:
                    consensus += f"[{AA_map[max_idx]}{AA_map[second_max_idx]}]"
            else:
                consensus += "."
                
    elif method.lower() == 'conservative':
        max_idx, max_value = np.argmax(profile.T, axis=1), np.amax(profile.T, axis=1)
        for idx, value in zip(max_idx, max_value):
            if value > cutoff:
                consensus += AA_map[idx]
            else:
                consensus += "."

    return consensus


def principal_component_analysis(features, labels, n_comp, location=None):

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
    if location is not None:
        fig.savefig(location + 'pca_features.eps', format='eps')
    plt.show()



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
        y_pred = model.predict(X[test])
        fpr, tpr, thresholds = roc_curve(
            y_true = y[test],
            y_score = y_pred,
            pos_label = 1
            )
#        viz = plot_roc_curve(model, X[test], y[test],
#                             name="", alpha=0, lw=3, ax=ax, color="royalblue")
        ax.plot(x = fpr, y = tpr, alpha=0, lw=3, color="royalblue")
        interp_tpr = np.interp(mean_fpr, fpr, tpr)
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
