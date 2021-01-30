---
layout: default
title: Methods
nav_order: 2
parent: Clustering
---

## Two Step

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

## MCL

