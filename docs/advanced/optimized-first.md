---
layout: default
title: Optimized clustering for large data sets
nav_order: 1
parent: Advanced
---

##  Optimizing the first clustering step

ClusTCR was designed with the aim of clustering large AIRR data sets exceeding millions of unique sequences. To accommodate for algorithmic runtime and memory limitations, we have introduced a batch clustering procedure.  Nonetheless, it is possible to achieve more performant results by changing the basic hyperparameters. This provides the user with even more flexibility, allowing them to decide upon the trade-off between clustering speed and accuracy.  

###  Adjusting supercluster size for speed improvement

The supercluster size (described by the `faiss_cluster_size` parameter) is an important determinant for the speed of ClusTCR. The MCL clustering step imposes a bottleneck to the speed of the algorithm. Therefore, by changing the `faiss_cluster_size` parameter, one can in- or decrease the total number of MCL processes. MCL is extremely efficient at clustering small groups of sequences, but becomes much slower when the size of the input data increases. In some cases, it may therefore be more beneficial to tune down the size of the superclusters. The `faiss_cluster_size` hyperparameter  can be set as a variable in the `Clustering` object. 

Let's illustrate this with an example:

```python
from clustcr import metarepertoire, Clustering

# Generate a sample of 1 million sequences
data = metarepertoire(directory='/emerson-2017/', data_format='immuneaccess')
```

We will now run ClusTCR with different preset values for the supercluster size.

```python
clustering_1 = Clustering(faiss_cluster_size=5000, n_cpus=16)
clustering_2 = Clustering(faiss_cluster_size=3000, n_cpus=16)
```

And observe the differences in algorithmic runtime.

```python
clustering_1.fit(data)
> 199.008 seconds # ~3min19sec

clustering_2.fit(data)
> 162.839 seconds # ~2min43sec
```

### Adjusting supercluster size for batch clustering

As discussed in the section on [batch clustering](../clustering/batch-clustering), the appropriate size for the training sample can be calculated as  `1000 * (fitting_data_size / faiss_cluster_size)`. When working with extremely large data sets, even the training sample may be too large to fit into memory. This problem can be accounted for by increasing the `faiss_cluster_size` parameter. Conversely, you can decrease the  `faiss_cluster_size`  to improve the speed of the batch clustering procedure. Suppose our dataset contains 10,000,000 unique CDR3 sequences, the optimal sample size is

<img src="https://latex.codecogs.com/svg.latex?\Large&space;S=1000\times\frac{10,000,000}{5,000}=2,000,000" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />

according to the default settings. By adjusting the  `faiss_cluster_size` we increase the total clustering time, but we won't need such a large training sample:

<img src="https://latex.codecogs.com/svg.latex?\Large&space;S=1000\times\frac{10,000,000}{10,000}=1,000,000" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />

