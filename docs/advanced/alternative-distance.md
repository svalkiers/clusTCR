---
layout: default
title: Alternative distance metric
nav_order: 2
parent: Advanced
---

##  Alternative distance metric

By default, distances between sequences are calculated by ClusTCR using the Hamming distance (HD) metric. The use of the HD implies that sequences within a single cluster all have identical length. However, there are known cases where CDR3 sequences of different length exert specificity towards the same epitope. Therefore, ClusTCR also offers the use of the Levenshtein distance (LD). To use LD instead of HD, you can change the default setting of the `distance_metric` hyperparameter to *levenshtein*:

```python
clustering = Clustering(distance_metric='levenshtein')
```

From our benchmarking, the use of LD over HD results in slightly higher cluster retention. For other clustering metrics (purity, fraction of clusters with >90% purity, consistency), no differences were observed.

**Note 1:** Calculating LD takes significantly longer than HD, especially for larger samples. We recommend only using this feature for data sets consisting of < 50,000 unique sequences.

**Note 2:** By using the LD, some clusters will contain sequences of unequal length. The `.summary()` and `.compute_features()` methods assume identical length of all the sequences within a cluster. Therefore, these methods should **not** be used when the distance metric is set to *levenshtein*.
