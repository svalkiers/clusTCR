---
layout: default
title: Methods
nav_order: 2
parent: Clustering
---

## Faiss

Implements a CDR3 clustering method using the [Faiss Clustering Library](https://github.com/facebookresearch/faiss).
By first converting the sequences into vectors using their chemical properties (more specifically, the mutation stability and z-scores), 
we can use the very efficient faiss K-means implementation to find a given amount of cluster centroids.
Using a faiss index, each vectorized sequence is assigned to a centroid, creating the clusters.

#### Characteristics 
- Very fast, especially when generating superclusters
- Scales up really well to extremely large datasets
- Average to low quality, when looking at purity and consistency
- Full retention, all sequences are clustered

#### Use

```python
faiss_clustering = Clustering(method='faiss')
output = faiss_clustering.fit(cdr3)
```

## MCL

MCL is a graph clustering algorithm that identifies dense network substructures in a graph or network by simulating stochastic flow. Hence, the algorithm performs a random walk on the graph, which is calculated using Markov chains. Essentially, MCL performs two operations on the stochastic matrix: expansion and inflation. Expansion is taking the power of the stochastic matrix. This determines how much flow is allowed between different regions of the graph. The inflation parameter is responsible for strengthening current between (already) strong neighbours, while at the same time weakening current between weak or distant neighbours.

MCL requires a graph as input. The nodes or vertices in this graph are represented by the CDR3 sequences, and the edges between them represent a maximum amino acid edit distance of one. More specifically, clusTCR uses the Hamming distance (HD), which implicitly assumes equal length of two sequences. A hashing function is used to bin sequences into hashes. This way, clusTCR only has to do pairwise comparisons of sequences within each individual hash. This drastically limits the total amount of pairwise comparisons, thus speeding up computation and scalability.

To learn more about this graph clustering algorithm, visit the [MCL homepage](https://micans.org/mcl/).

#### Characteristics 

- Fast for smaller datasets (< 50K sequences)
- Poor scales up to large datasets
- High quality, when looking at purity and consistency
- Maximum edit distance implies incomplete retention. Only a part of the sequences are clustered.

#### Use

```python
mcl_clustering = Clustering(method='mcl')
output = mcl_clustering.fit(cdr3)
```

## Two Step

Two-step clustering algorithm combines the speed of the Faiss, combined with the accuracy of the Markov clustering algorithm (MCL). During the first clustering step we use Faiss' efficient K-means implementation to rapidly subdivide data sets of CDR3 sequences into superclusters. In the next step, MCL is applied on each individual supercluster to identify groups of epitope-specific CDR3 sequences. **Two-step is clusTCR's default clustering method.**

#### Characteristics 

- Very fast
- Scales up very well to large datasets
- High quality, when looking at purity and consistency
- Maximum edit distance of MCL step implies incomplete retention. Only a part of the sequences are clustered.

#### Use

```python
ts_clustering = Clustering(method='two-step')
output = ts_clustering.fit(cdr3)
```