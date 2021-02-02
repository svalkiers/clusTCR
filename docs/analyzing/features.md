---
layout: default
title: Features
nav_order: 1
parent: Analyzing
---


## Features

Clusters can be represented by a feature matrix that describes several properties of the amino acid sequence within that cluster, including a number of physicochemical properties, entropy, size, length and generation probability. **Table 1** contains a list of all features calculated by clusTCR, and provides a description for all of them.

**Table 1: Computed cluster features**

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

These features can be calculated by creating a `FeatureGenerator` object. The input provided should correspond to the output of the clustering procedure. From this object, you can compute all features provided in **Table 1**.

```python
from clustcr import FeatureGenerator
fg = FeatureGenerator(clusters)
features = fg.compute_features()
```

**NOTE:** calculating generation probabilities is time-consuming. You will receive a warning asking you whether you want to proceed.

```
Calculating generation probabilities may take a while. Are you sure you want to continue?

Confirm: [Y/N] 
```

If you select `N` the algorithm will generate a feature matrix nevertheless. However, keep in mind that this matrix will not contain values for *pgen_avg* and *pgen_var*.