---
layout: default 
title: Batch Clustering
nav_order: 1 
parent: Large Datasets
grand_parent: Clustering
---



## Batch Clustering Large Datasets

Is your CDR3 dataset so big it simply cannot fit into your RAM memory when loading it? Or do you have a large data set, consisting of a large number of samples? You've come to the right place.

ClusTCR includes methods to cluster data in batches with a comparable accuracy to the normal way of clustering. In
this tutorial we'll go over the steps.

#### Prepping the data

The batch clustering functionality provided within the ClusTCR package requires that batches are separated into different files. Thus, if the data you want to cluster is stored in one file, you should split it up into different files with batch sizes that comfortably fit into your device's memory. Here we give an example of how you could do this. Suppose you have a very large repertoire file which you would like to split up into different files, each containing 100,000 sequences. Here's how you do it. Using the `pandas` library, you can read in a large data frame in different chunks. We can use this functionality to write each chunk to an individual file. An example:

```python
import pandas as pd
import os

# Make a new directory to store the chunks
target_dir_name = 'my_tcr_chunks/'
os.mkdir(target_dir_name)
# Make a TextFileReader object to iteratively read the large repertoire file
chunk_reader = pd.read_csv('large_file.tsv', sep='\t', chunksize=100000)
# Use a for loop to write the different chunks to distinct files
for n, chunk in enumerate(chunk_reader):
    chunk.to_csv(os.path.join(target_dir_name, 'my_tcr_chunk_%s' % n))
```

#### Training

The first part of the Faiss clustering consists of training. In this step the centroids of the clustering
method are calculated, to which we later assign all the sequences and retrieve the superclusters.

3 parameters are necessary for training with a subset of the fitting data.

| parameter | explanation |
|:-------------|:------------------|
| faiss_training_data | the sample used for training |
| fitting_data_size | the amount of CDR3 sequences that will be clustered |
| max_sequence_size | the size of the largest sequence (the faiss method works on vectorized sequences and naturally expects these vectors to be of the same size) |

Normally, the whole dataset is used for the calculation of the centroids. In a case of restricted memory, this is
unfortunately not always possible. What we do instead is take a meaningful sample of our dataset for training. To
achieve comparable accuracy to the standard way, the size of this sample should be
around `1000 * (fitting_data_size / faiss_cluster_size)` where `faiss_cluster_size` is also a parameter of the `Clustering` interface as described [here](how-to-use). 
Any larger and there won't be a noticeable difference.


```python
from clustcr import Clustering
clustering = Clustering(faiss_training_data=sample,
                        fitting_data_size=fitting_data_size,
                        max_sequence_size=max_sequence_size)
```


#### Batch Preclustering

Next up, we do our "preclustering".
In other words, we assign all the sequences to a centroid.
Because our data is too large to load at once, we do this preclustering in batches.

```python
for file in files:
    # Load your data
    data = process_file(file)
    clustering.batch_precluster(data)
```

To store the preclusters, ClusTCR stores them on disk, with a file for each precluster.


#### Batch Clustering

Lastly, the MCL clustering is performed on each of the preclusters.
To make sure that memory restrictions are met, this is done in batches as well.
Using a generator function, data is loaded and clustered batch per batch.
To optimize this process further, a couple of preclusters are loaded at once and clustered using multiprocessing.


```python
for clusters in clustering.batch_cluster():
    # Do something with the clusters
    print(clusters.clusters_df)
```

At each iteration a ClusteringResult is returned.
To store the full result of the clustering it's recommended to store them on disk.

#### Batch Cleanup

In the process, a temporary directory is created that includes a file for each precluster.
To remove this at the end, the cleanup function can be called.

```python
clustering.batch_cleanup()
```

#### Applications

In the next section, we discuss how we can use the functionality of the batch clustering procedure to construct a cluster matrix that describes how the clusters are distributed across samples.

