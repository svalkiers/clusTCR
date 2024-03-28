---
layout: default
title: Features
nav_order: 1
parent: Analyzing
---


## Features

Clusters can be represented by a feature matrix that describes several properties of the amino acid sequence within that cluster, including a number of physicochemical properties, entropy, size, length and generation probability. **Table 1** contains a list of all features calculated by *clusTCR*, and provides a description for all of them.

**Table 1: Cluster features calculated by *clusTCR*, and their description.**

| feature               | description                                                  |
| :-------------------- | ------------------------------------------------------------ |
| h                     | Cluster entropy. This value describes the average information content per amino acid position (ignoring C and F at the first and last position respectively). as determined by the Shannon entropy at each position. A correction is applied that normalizes for cluster size. |
| size                  | Cluster size. Number of sequences in the cluster.            |
| length                | Length of CDR3 sequences in the cluster. Due to the fact that we use Hamming distance, (HD) all sequences within a cluster have exactly equal lengths (HD assumes equal length of sequences). |
| basicity_avg          | Average basicity of sequences in the cluster.                |
| basicity_var          | Varience in basicity of sequences in the cluster.            |
| hydrophobicity_avg    | Average hydrophobicity of sequences in the cluster.          |
| hydrophobicity_var    | Varience in hydrophobicity of sequences in the cluster.      |
| helicity_avg          | Average helicity of sequences in the cluster.                |
| helicity_var          | Varience in helicity of sequences in the cluster.            |
| mutationstability_avg | Average mutation stability of sequences in the cluster.      |
| mutationstability_var | Varience in mutation stability of sequences in the cluster.  |
| pgen_avg              | Average generation probability of CDR3 sequences in the cluster. Generation probability is calculated using the *olga* module. |
| pgen_var              | Variance in generation probability within the cluster.       |

**Note:** if you want to compute the generation probability for the TCR alpha chain, user must specify the `chain='A'` parameter in the `Clustering` object.

These features can be calculated by calling the `compute_features()` function on a `ClusteringResult` object (see [clustering section](../clustering/how-to-use)). The code block below shows a brief example of a workflow for calculating cluster features.

```python
from clustcr import datasets, Clustering

# Load some data
data = datasets.test_cdr3()

# Perform clustering
clustering = Clustering(chain='B') # change to chain='A' for alpha chain
output = clustering.fit(data)

# Compute features of ClusteringResult object
features = output.compute_features(compute_pgen=True)
```

Calculating generation probabilities (*pgen*) is time-consuming. A `compute_pgen` parameter is provided, so the user can indicate whether they want to calculate *pgen* values. If this parameter is set to `False`, a feature matrix *without* pgen values is computed. Note that *pgen* is necessary if you want to use the cluster quality classifier's functionality (see section on [exploring clustering results](exploration)).