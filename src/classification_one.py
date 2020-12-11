# -*- coding: utf-8 -*-
"""
author: Sebastiaan Valkiers
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from networkTCR import Dataloader, Clustering, Features, Metrics
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor
from sklearn.decomposition import PCA
from sklearn import metrics

# Load data
Data = Dataloader("../data/vdjdb_trb.tsv")
vdjdb = Data.read_bm_file(q=0)

# Perform clustering
Clust = Clustering(set(vdjdb["CDR3"]))
edges = Clust.create_network()
nodes = Clust.network_clustering(edges)

# Calculate features of clusters
Feat = Features(nodes)
aavar = Feat.calc_variation()
pchem = Feat.calc_physchem()
pgen = Feat.calc_pgen()

features = Feat.combine(aavar, pchem, pgen) # combines results into joint DataFrame

# Generate labels (purity)
Metr = Metrics(nodes, vdjdb)
cm = Metr.calc_confmat()[0]
labels = []
for i in cm:
    labels.append(cm[i].max() / np.sum(cm[i]))
    
features["label"] = labels
features.dropna(inplace=True)

X = np.array(features.iloc[:,:-1])
y = np.array(features.iloc[:,-1])

X = StandardScaler().fit_transform(X)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.2)

regr = RandomForestRegressor(n_estimators=100,
                              n_jobs=-1)
regr.fit(X_train, y_train)

pred = regr.predict(X_test)

print('Mean Absolute Error (MAE):', metrics.mean_absolute_error(y_test, pred))
print('Mean Squared Error (MSE):', metrics.mean_squared_error(y_test, pred))
print('Root Mean Squared Error (RMSE):', np.sqrt(metrics.mean_squared_error(y_test, pred)))
mape = np.mean(np.abs((y_test - pred) / np.abs(y_test)))
print('Mean Absolute Percentage Error (MAPE):', round(mape * 100, 2))
print('Accuracy:', round(100*(1 - mape), 2))

fig, ax = plt.subplots(figsize=(12,8))
df = pd.DataFrame({"feature":list(features.columns[:-1]), 
                   "importance":regr.feature_importances_})
df.sort_values(by="importance", ascending=False, inplace=True)
sns.barplot(x=df["importance"],
            y=df["feature"], 
            orient="h",
            color=sns.color_palette()[0])
ax.set_xticklabels(fontsize=16, labels=np.round(np.arange(0,0.26,0.05),2))
ax.set_yticklabels(fontsize=16, labels=df["feature"])
ax.set_xlabel("Importance", fontsize=20)
ax.set_ylabel("")
ax.set_title("Feature importance - VDJdb", fontsize=26)

pca = PCA(n_components=5)
x_new = pca.fit_transform(X)

fig, ax = plt.subplots(figsize=(12,8))
score = x_new[:,0:2]
xs = score[:,0]
ys = score[:,1]
coeff = np.transpose(pca.components_[0:2, :])
n = coeff.shape[0]
labels = features.columns[:-1]
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


