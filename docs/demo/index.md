---
layout: default
title: Demo
nav_order: 8
has_children: False
---

## clusTCR demo

In this section we will go over some typical scenarios that you might encounter when you will use clusTCR.

### > Scenario 1: clustering a single repertoire file

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

```python
print(data)
>> 0        CATTGTSGGPSQSTQYF
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

|       | CDR3             | cluster |
| :---- | :--------------- | :------ |
| 0     | CASSYSRASSGNTIYF | 0       |
| 1     | CASSYSRQSSGNTIYF | 0       |
| 2     | CASSYSGQGSGNTIYF | 0       |
| 3     | CASSYSRQGSGNTIYF | 0       |
| 4     | CASSYSRPSSGNTIYF | 0       |
| ...   | ...              | ...     |
| 25065 | CASSYGSSSTDTQYF  | 4516    |
| 25066 | CASSPQWATGNTIYF  | 4517    |
| 25067 | CSSSPQWATGNTIYF  | 4517    |
| 25068 | CATSRDQGGYNEQFF  | 4518    |
| 25069 | CATSRDRGGYNEQFF  | 4518    |

...or directly save the clusters to a file using the `write_to_csv()` method. This will require you to specify a path to where the results should be saved, else they will be stored in the current working directory with the generic *clusTCR_clusters.csv* file name.

```python
result.write_to_csv(path='/results_folder/myclusters.csv')
```

Alternatively, you can just get a flat notation of the clustering result in the form of a list, in which each element is a list of the sequences within a cluster.

```python
result.cluster_contents()
>> [
	['CASSRWTGTNTGELFF', 'CASSSWTGTNTGELFF'],
	['CASSLPGQGMNTEAFF', 'CASSLVGQGMNTEAFF', 'CASSLVGLGMNTEAFF'], ...
]
```

If you want a more concise overview of the clustering results, you can call the `summary()` method on the `ClusteringResult` object. 

```python
result.summary()
```

|      | cluster_idx | size | motif               |
| :--- | :---------- | :--- | ------------------- |
| 0    | 2238        | 245  | CASS[RS]WTGTNTGELFF |
| 1    | 2249        | 241  | CASSLVGQGMNTEAFF    |
| 2    | 1252        | 237  | CASSLALQ[RG]YGNTIYF |
| 3    | 3329        | 231  | C[AS]SSGARLGYREKLFF |
| 4    | 1246        | 226  | CASSYS[RK]GGAGIWAFF |
| ...  | ...         | ...  | ...                 |
| 4514 | 1593        | 2    | CASS[FT]TTGGGNEQFF  |
| 4515 | 3640        | 2    | CASSPPRG[RQ]GETQYF  |
| 4516 | 1577        | 2    | CASSSYDRK[AV]YEQYF  |
| 4517 | 3624        | 2    | CASSFSGT[GL]GNTIYF  |
| 4518 | 0           | 2    | CASSF[GS]GGAGDEQFF  |



### > Scenario 2: clustering a set of repertoires simultaneously

Suppose that you have a data set containing various TCR repertoire samples, which you want to cluster simultaneously. Most likely, the size of this data will be larger than your maximum RAM capacity. To solve this problem, we can use clusTCR's batch clustering functionality.

Start by importing the different modules that we'll need:

```python
from clustcr import read_cdr3, metarepetoire, Clustering
import os
```

First, make sure your files are stored in a separate directory. Batch clustering requires three specific parameters:

- `faiss_training_data`: A training sample, from which it will compute the cluster centroids.
- `fitting_data_size`: The total number of sequences in the data set.
- `max_sequence_size`: The length of the largest CDR3 sequence in the data set.

We will start by determining the size of the data set. Then, we can calculate the recommended training sample size. For the sake of this example, we will illustrate this procedure for an immuneACCESS data set.

```python
# First, we define the path to the data directory
datadir = '/path_to_data/'

# Now we count the number of sequences that are present in the data set
# We use os.listdir() to list all the files in the specified directory
total_sequences = 0
for file in os.listdir(datadir):
    total_cdr3s += len(read_cdr3(datadir + file, 
                                 data_format='immuneaccess'))
```

From this number, we can calculated the recommended sample size (see [clustering large data](../clustering/large-data.md) section). We can then take a sample from the data set using the `metarepertoire()` function.

```python
training_sample_size = round(1000 * (total_sequences / 5000))
training_sample = metarepertoire(directory=datadir,
                                 data_format='immuneaccess',
                                 n_sequences=training_sample_size)
```

We assume that the largest sequence in the sample also corresponds to the maximum sequence size in the complete data set. If we allow this assumption, the maximum sequence size can be easily determined:

```python
max_seq_len = training_sample.str.len().max()
```

Next, we make a Clustering object. However, this time we need to specify the batch clustering-specific parameters.

```python
clustering = Clustering(faiss_training_data=training_sample,
                        fitting_data_size=total_sequences,
                        max_sequence_size=max_seq_len,
                        n_cpus=8) # Multiprocessing using 8 CPUs
```

Finally, we cluster the sequences, starting first with the batch pre-clustering...

```python
for file in files:
    # Load your data
    data = read_cdr3(file=os.path.join(datadir, file),
                    data_format='immuneaccess')
    clustering.batch_precluster(data)
```

...followed by the actual clustering step:

```python
for cluster in clustering.batch_cluster():
	# Do something with the clusters
	# For example: - print(cluster.clusters_df)
	#			   - cluster.write_to_csv()
	#			   - ...
```

Once finished, you can clean up the intermediate results using the following command:

```python
clustering.batch_cleanup()
```

