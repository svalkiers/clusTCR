# TCR Clustering

This library aims to do a fast clustering of T-cell receptors using a combination of chemical properties as features.  
It relies on the [faiss](https://github.com/facebookresearch/faiss/) library for clustering.

This project is made during an Internship at the [Adrem Data Lab](https://www.uantwerpen.be/en/research-groups/adrem-data-lab/) 
of the University of Antwerp with supervision and help of [Pieter Meysman](https://www.uantwerpen.be/nl/personeel/pieter-meysman/) and others.

 
## Install
 
 The first step to be able to run this project is installing the faiss library.
 For this they have an [instruction manual](https://github.com/facebookresearch/faiss/blob/master/INSTALL.md).
 The easiest way to install is using conda which comes down to
 
 ```
 conda install faiss-cpu -c pytorch
 ```
 
 If the installation fails because of an incorrect python version, make a new conda python environment with a version that is supported.
 This should fix the problem.
 
 

## Clustering

A simple clustering of CDR3 sequences can be made by supplying a pandas.Series with the sequences

```python
from faiss_clustering import FaissClustering
cdr3 = data['CDR3']
clustering = FaissClustering.cluster(data)
```

### Clustering Output
To analyse the resulting clustering, we can look at the clusters
 ```python
clustering = FaissClustering.cluster(data)
for cluster in clustering.get_cluster_contents():
    print(cluster)  # Prints for example ['CAT', 'CAI', 'CTT']
``` 

### Properties

The above clustering is done using an optimal combination of chemical properties.
However, the properties can be manually given:
```python
from faiss_clustering import properties
clustering = FaissClustering.cluster(data, properties=[properties.BASICITY, properties.ISOELECTRIC])
``` 

### Items per cluster

The method generates clusters of an average size of 10, this can also be adjusted:
 ```python
clustering = FaissClustering.cluster(data, items_per_cluster=30)
``` 


## Distance Calculation

Calculating the edit distance for each possible pair in large datasets quickly becomes extremely time consuming.
By first clustering the sequences, we can calculate the distance only for pairs in the same cluster to speed up the process
while still maintaining decent accuracy.

```python
from faiss_clustering import DistancePairs
clustering = FaissClustering.cluster(data)
distance_pairs = DistancePairs.generate(clustering, data)
```

### Distance Output

From the distance pairs we can easily output a dictionary where pairs are grouped per distance

```python
distance_pairs.to_dict()
# For example { 1: [('CAI', 'CAT'),..], 2: [..] }
```

### Accuracy

Because we only calculate the distance for the pairs in a cluster, we won't be able to find all the pairs that are close by.
To check how well we do, we compare the amount of pairs with distance 1 we have found against the total amount of pairs with distance 1, for 
which there is a fast hashing algorithm.

```python
distance_pairs.percentage_of_distance1_pairs_found()
```


## Input

To process a dataset, we can import it with ease using the `process_csv` function.
We simply provide the files and the cdr3 column name, and it will import it, drop duplicates and drop incomplete sequences.

```python
from faiss_clustering import process_csv
cdr3 = process_csv('filename.csv', 'cdr3')
# cdr3 is a pandas.Series
```

The same function can be used to process multiple files, where they are added together

```python
cdr3 = process_csv(['filename1.csv', 'filename2.csv'], 'cdr3')
```

### Datasets with epitopes

When epitopes are also present in a dataset, we simply provide the column name and the function will return the epitopes as a series as well

```python
cdr3, epitopes = process_csv('filename.csv', 'cdr3', epitope_column_name='Epitope')
```

