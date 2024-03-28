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
| chain | *A* or *B*<br />Specify alpha (*A*) or beta (*B*) chain. This choice does not influence the clustering process, but the information will be used during the downstream cluster analysis. Default = *B*. |  |
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

#### Including CDR3 alpha chain
In case you have data that includes the alpha chain, you can optionally use this in the clustering (beware that this will drastically change your results).
To include, simply pass a Series of your alpha chain to the fit method.
As an example, we provide a vdjdb dataset that includes this data.

```python
data = datasets.vdjdb_paired()
cdr3, alpha = data['CDR3_beta'], data['CDR3_alpha']
output = clustering.fit(cdr3, alpha=alpha)
```

 #### [BETA] Including V gene information

In addition to CDR3-based clustering, ClusTCR provides the ability to include V gene information into the clustering process. When V gene clustering is enables, TCR sequences will first be sorted by V gene family, and clustering will be applied within each group of sequences that belong to that V gene family. By doing so, clustering accuracy will be increased at the cost of clustering retention (i.e. less sequences will end up in a cluster).

You can include V gene information into the clustering process by setting the `include_vgene = True`  in the `.fit()` method. In addition, this process requires the user to specify the names of the columns containing the CDR3 and V gene information. Below, you can find an example of how this works in practice:

```python
import pandas as pd
from clustcr import Clustering, datasets

# Import a file that contains at least a CDR3 column and a V gene column
data = pd.read_csv("mytcrfile.csv")
# Initiate a Clustering object
clustering = Clustering()
# Include V gene information by setting the include_vgene parameter to True
output = clustering.fit(
    data,
    include_vgene = True, # Enable V gene clustering
    cdr3_col = "cdr3", # Specificy CDR3 column name
    v_gene_col = "vgene" # Specificy V gene column name
	)
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

You can explore the clustering results by executing the `.summary()` on the ClusteringResult object.  This will provide you with a `pandas.DataFrame` that contains the cluster index, number of sequences in the cluster, and a consensus motif for that cluster.  For every position, the amino acid frequency is calculated. The most dominant amino acid is selected and if it exceeds the predefined **cut-off** (default = 0.7), that position is represented by the dominant amino acid in **upper case**. Else, if the sum of the frequencies of the two most dominant amino acids exceed the cut-off, both of them are considered. However, if the frequency of one the two is 2x larger than the other, the most dominant amino acid of the two will be used to represent that position. To emphasize that the frequency of this amino acid on itself does not exceed the cut-off, it will be shown in **lower case**. Otherwise shared dominance of amino acids is indicated by **square brackets (' [ ] ')**. Positions where neither criteria are met are indicated with a **wild card symbol (' . ')**.

```python
output.summary()
```

|   |  cluster_idx | size  |              motif|
| :------- | :------- | :---------- | : ------- |
| 0 |           25 |  24 | CASSgg.YGYTF        |
|1  |           0   | 23  |CASS.RSTDTQYF|
|2  |          15  |  15  | CASSEA[AS]GGFYNEQFF |
|3  |          41  |  15  | CASSL[LM]GPGQPQHF |
|4  |           5  |  14  | CSAR.GLNNEQFF |
|.. |          ...  | ...  |                ...|
|235 |        135  |   2  | CASS[LP]GWGLDQPQHF |
|236 |        137  |   2  | CASSLLGQ[DY]NSPLHF |
|237 |         26  |   2  | CASSLEG[DY]TEAFF |
|238 |        139  |   2  | CASSSTGGGG[AT]EAFF  |
|239 |       239  |   2  |CASS[EL]GRETQYF|