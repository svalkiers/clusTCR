---
layout: default 
title: Feature Matrix
nav_order: 2
parent: Large Datasets
grand_parent: Clustering
---

## Multi Repertoire Feature Matrix

When combining TCR repertoires from multiple people into a metarepertoire, we lose sight of each person's unique data.
To solve this, we generate a feature matrix to be able to compare data from different sources.

Here's an example matrix

|  | 1 | 2 | 3 | 4 | ... | 
|:--- |:--- | :--- | :--- | :--- | :--- | 
| HIP00707.tsv | 29 | 38 | 16 | 8 | ... |
| HIP00710.tsv | 16 | 8 | 17 | 4 | ... | 
| HIP13923.tsv | 2 | 3 | 0 | 1 | ... |
| HIP14106.tsv | 0 | 2 | 0 | 0 | ... | 

Each row in the matrix is one person's data and each column is a cluster (corresponding to the cluster-id from the clustering output).
The values in the matrix are the amount of TCRs of a person found in a cluster.
For example, there are 29 TCRs from *HIP00707.tsv* found in cluster 1.


### Implementation

To achieve this matrix, we use the [batch clustering](batch-clustering) methods with a couple of added tweaks.
After initializing the *Clustering* object like normal, we perform the batch preclustering.
Now however, we're adding a name for every dataset we precluster. This name will be used in the feature matrix.

In our case, we simply used the filenames

```python
for filename, data in loaded_datasets:
    clustering.batch_precluster(data, name=filename)
```

Afterwards, we do the batch clustering and simply specify that we want the feature matrix to be calculated

```python
for cluster in clustering.batch_cluster(calc_feature_matrix=True):
    print(cluster.clusters_df)
```

Lastly, the feature matrix can be accessed

```python
clustering.batch_feature_matrix()
```








