---
layout: default
title: Clustering
nav_order: 3
has_children: true
---


## Two Step Clustering


To perform clustering, import the `Clustering` module from *clusTCR*. A test set of CDR3 sequences is also provided within *clusTCR*, which will be used in this tutorial.

```python
from clustcr import clusTCR
cdr3 = clusTCR.test_cdr3()
```

The next step is to generate a `clusTCR.Clustering` object. The `Clustering` method takes three arguments:

1. A `pandas.Series` of CDR3 amino acid sequences.
2. The desired clustering method. Currently, *clusTCR* contains three clustering approaches: **MCL-based clustering, faiss-based clustering, and a two-step clustering** procedure which combines the first two methods for an increase in speed and accuracy. Valid options are: *mcl*, *faiss* or *two-step*. We recommend using *mcl* for data sets containing < 50,000 CDR3 sequences, and *two-step* for all data sets with > 50,000 sequences.
3. Number of CPUs. This drastically increases the speed of *clusTCR*. By default, this value is set to 1, but we recommend using at least 4 CPUs. If you want to get the maximum out of your machine, but you don't known the total amount of cores, you can use `n_cpus = 'max'`.

To perform the actual clustering, you will need to execute the `get_clusters()` method on the newly created object.

```python
from clustcr import Clustering

clustering = Clustering(cdr3, method='two-step', n_cpus=4)
output = clustering.get_clusters()
```
