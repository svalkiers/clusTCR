---
layout: default
title: Usage
nav_order: 1
parent: Clustering
---


## Clustering

Our interface follows the `sklearn` practices, so it might feel familiar.  
To perform clustering, import the `Clustering` class from *clustcr*.

```python
from clustcr import Clustering
clustering = Clustering()
```

Any parameters should be passed to the `Clustering` object at creation, the following are available

| parameter | explanation | default |
|:-------------|:------------------|:------|
| method |  *mcl*, *faiss* or *two-step*. <br> We recommend using *mcl* for data sets containing < 50,000 CDR3 sequences, and *two-step* for all data sets with > 50,000 sequences. For more information check out the [methods page](methods). | two-step  |
| n_cpus | Number of CPUs used in the MCL clustering. This drastically increases the speed of *clusTCR*. When set to 'all', all of your CPUs will be used. | 1  |
| use_gpu | Usage of GPU in the Faiss Clustering training step (needs clustcr-gpu to be installed) | False |
| faiss_cluster_size | The size of the clusters that faiss will generate, either using the faiss or the two-step method. | 5000 |
| mcl_params | MCL hyperparameters, which should be a list of \[inflation, expansion\] | \[1.2, 2\]  |
| faiss_training_data, fitting_data_size, max_sequence_size | Only used for clustering in batches, see [clustering large data](large-data) | /  |

To perform the clustering, our interface provides a `fit` method which expects a `pandas.Series`.
This method can be called multiple times, each time returning the clustering result.  
For this tutorial, we use the `test_cdr3` dataset.

```python
from clustcr import datasets
cdr3 = datasets.test_cdr3()
output = clustering.fit(cdr3)
```

### ClusteringResult

#### Dataframe

A dataframe containing the clusters can be accessed 

```python
output.clusters_df
```

|      |           CDR3  | cluster |
| :------- | :------- | :---------- |
| 0   |   CASSPSGTPYEQYF |       0 |
| 1   |   CASSPSGTPYERYF |       0 | 
| 2   |  CASNELASGTDTQYF |       1 |
| 3   |  CASSELASGTDTQYF |      1 |
| 4   |  CASSALASGTDTQYF |       1 |
| ..  |              ... |     ... |
| 637 |  CASSPRTSGTYEQYF |    199 |
| 638 | CASSFTLGTGGVEQYF |     200 |
| 639 | CASSITLGTGGVEQYF |     200 |
| 640 | CASSLIGVSSYNEQFF |     201 |
| 641 |  CASSLRGVSSYNEQFF |      201 |

#### CSV

To quickly store the clusters to file, the `write_to_csv` method can be used.
A path is optional, by default *clusTCR* will save it in the current directory.
```python
output.write_to_csv()
```

#### Cluster Contents

To have a representation of the CDR3s in each cluster, the following method can be used
```python
output.cluster_contents()
>> [
    ['CASSPSGTPYEQYF', 'CASSPSGTPYERYF'], 
    ['CASSFTLGTGGVEQYF', 'CASSITLGTGGVEQYF'], ...
]
```



#### Summary

You can explore the clustering results by executing the `.summary()` on the ClusteringResult object. 
This will provide you with a dataframe that contains the cluster index, number of sequences in the cluster, and a consensus motif for that cluster. 
Shared dominance of amino acids is indicated by square brackets (' [ ] ').
If there is no dominant amino acid, the position is denoted by a single dot (' . ').

```python
output.summary()
```

|   |  cluster_idx | size  |              motif|
| :------- | :------- | :---------- | : ------- |
| 0 |            50 |   17 |   CASSPSGTPYE[RQ]YF|
|1  |           5   | 14   |CASS[EA]LASGTDTQYF|
|2  |           17  |  13  |      CASTRGSTDTQYF|
|3  |           35  |  13  |       CASSSANYGYTF|
|4  |           18  |  12  |       CASSSLNTEAFF|
|.. |          ...  | ...  |                ...|
|197 |          86  |   2  |     CASTGG[LS]GYTF|
|198 |          85  |   2  |   CASS[LF]GSYNEQFF|
|199 |          84  |   2  | CASSPRTSG[PT]YEQYF|
|200 |          83  |   2  |CASS[IF]TLGTGGVEQYF|
|201 |         201  |   2  |CASSL[RI]GVSSYNEQFF|
