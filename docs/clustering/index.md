---
layout: default
title: Clustering
nav_order: 3
has_children: true
---


## Two Step Clustering


To perform clustering, import the `Clustering` module from *clustcr*. A test set of CDR3 sequences is also provided within *datasets*, which will be used in this tutorial.

```python
from clustcr import datasets
cdr3 = datasets.test_cdr3()
```

The next step is to generate a `clustcr.Clustering` object. The `Clustering` object takes in all possible options:

1. "method": The desired clustering method. Currently, *clusTCR* contains three clustering approaches: **MCL-based clustering, faiss-based clustering, and a two-step clustering** procedure which combines the first two methods for an increase in speed and accuracy. Valid options are: *mcl*, *faiss* or *two-step*. We recommend using *mcl* for data sets containing < 50,000 CDR3 sequences, and *two-step* for all data sets with > 50,000 sequences.
2. "n_cpus": Number of CPUs. This drastically increases the speed of *clusTCR*. By default, this value is set to -1, which means all of your CPUs will be used.
3. "faiss_cluster_size": The size of the clusters that faiss will generate, either with the faiss method or the two-step.
4. "mcl_params": MCL Hyperparameters \[inflation, expansion\]

To perform the actual clustering, you will need to execute the `fit()` method on the newly created object, which takes in a `pandas.Series` of CDR3 amino acid sequences..


```python
from clustcr import Clustering

clustering = Clustering(method='two-step', n_cpus=4)
output = clustering.fit(cdr3)
```

The resulting output is a ClusteringResult object, where a dataframe containing the clusters can be accessed 

```python
print(output.clusters_df)
```

