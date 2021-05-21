---
layout: default
title: Exploration
nav_order: 2
parent: Analyzing
---


## Explore clustering results

### Initiate the ClusterAnalysis object

clusTCR provides a number of options for exploration of clustering results. To perform these analysis, you must initiate a `ClusterAnalysis` object, which takes cluster features as its only argument.

```python
from clustcr import ClusterAnalysis
analysis = ClusterAnalysis(features)
```

#### PCA

To get a quick overview of the newly generated features, you can perform a principal component analysis (PCA). clusTCR contains a built-in PCA functionality, which can be conveniently executed by performing the `._pca()` method on a `ClusterAnalysis` object.

```python
analysis.pca(features)
```

Performing a PCA with clusTCR will provide a figure of the PCA loadings.

<p align="center">
    <img src="cluster_features_PCA.png" alt="drawing" width="500"/>
</p>




#### Analyzing cluster quality

An additional feature of clusTCR's analysis module is predicting the quality of a cluster. Here, quality is determined as the **purity** of a cluster. This feature is particularly useful when no information is available about the target epitopes of the clustered TCR sequences. We trained a classification model that predicts whether an individual cluster will be of good (1) or bad (0) quality. Good clusters have a predicted purity of > 0.90, low quality clusters have a purity < 0.90.

##### Using the pre-trained model

```python
predictions = analysis.predict_quality()
```

##### Training your own model

You can also train your own model by creating a `ModelTraining` object. To do this, you will need to provide the following three arguments: features (see [features](./features.md) page), results from the clustering and epitope data corresponding to the CDR3 sequences used for clustering.

```python
from clustcr import ModelTraining
model = ModelTraining(features, clusters, epitopes)
clf = model.fit()
```

You can evaluate your own model using the `._evaluate()` method. This will perform 10-fold stratified cross-validation and outputs a receiver operating characteristic (ROC) curve with corresponding area under the curve (AUC) value.

```python
# Evaluate your own model
model.evaluate()
```

<p align="center">
    <img src="cluster_quality_ROC.png" alt="drawing" width="500"/>
</p>

Save your model using the following method:

```python
model.save(clf, '/path_to_model/my_custom_model.pkl')
```

If you want to use this model to make quality predictions, just specify the path to your model in the *model* parameter of the `_predict_quality()` method.