---
layout: default
title: Demo
nav_order: 7
has_children: False
---

## clusTCR demo

In this section we will go over some typical scenarios that you might encounter when you will use clusTCR.

### Scenario 1: clustering a single repertoire file

Say you want to cluster an individual RepSeq file, you can simply load it in and perform clustering. Make sure you are using one of the supported input formats. For more info about input formats supported by clusTCR, see the [importing](../importing/index.md) section. Let's start by importing the appropriate modules:

```python
from clustcr import read_cdr3, Clustering
```

Next, load in your data. In this demo we will use the example of an immuneACCESS file. First, adjust the *file* parameter with the path to the file you want to cluster.

```python
data = read_cdr3(file='/path_to_data/HIP00110.tsv',
                data_format='immuneaccess')
```

This should return a `pandas.Series` of CDR3 sequences. Example:

```
0        CATTGTSGGPSQSTQYF
1           CASSLRVGGYGYTF
2          CASRRLGGLNTEAFF
3           CASSLRGSGNTIYF
4            CASRQDGSTEAFF
               ...        
77475            RVGSCEQYF
77476    CASSRYDLPGPRDTQYF
77477    CASSLVFYGQGQETQYF
77478           CATMGHGYTF
77479          CATSFRGEAFF
Length: 77480, dtype: object
```

Next, we will perform the actual clustering. First, we must define the clustering parameters (for more info, see the [Clustering](../clustering/index.md) section), creating a `Clustering` object. Then we fit our data onto the newly created `Clustering` object.

```python
clustering = Clustering(n_cpus=8) # Clustering parameters
result = clustering.fit(data) # This will generate the ClusteringResult object
```

That's it! You have successfully clustered your RepSeq file with clusTCR. To retrieve the clusters, you can call the `clusters_df` method on the `ClusteringResult` object...

```python
result.clusters_df
```



...or directly save the clusters to a file using the `write_to_csv()` method. This will require you to specify a path to where the results should be saved, else they will be stored in the current working directory with the generic *clusTCR_clusters.csv* file name.

```python
result.write_to_csv(path='/results_folder/myclusters.csv')
```

Alternatively, if you want a more concise overview of the clustering results, you can call the `summary()` method on the `ClusteringResult` object. 

```python
result.summary()
```

