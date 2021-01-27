---
layout: default
title: Evaluating
nav_order: 4
---


## Evaluating Clustering Quality

If you have epitope data available, you can evaluate the quality of the clustering output using the `Metrics` module. Therefore, start by importing epitope data corresponding to the clustered CDR3 sequences. *clusTCR* contains a data set of epitope targets corresponding to the sequences used in the previous section. The epitope data should always have the following structure: two tab-separated columns with headers `['CDR3', 'Epitope']`.

```python
from clusTCR import test_epitope
epitope = test_epitope() # Epitope data corresponding to CDR3 sequences
```

Currently, *clusTCR* provides the following metrics for cluster evaluation:

- **Retention**: Fraction of sequences that is assigned to *any* cluster.
- **Purity**: Fraction of sequences within a cluster targeting the same epitope.
- **Consistency**: Fraction of sequences targeting the same epitope that are assigned to the same cluster.

```python
from clusTCR import Metrics

metrics = Metrics(output, epitope) # Provide clustering output and epitope data
retention = metrics.retention()
purity = metrics.purity()
consistency = metrics.consistency()
```

Printing these variables will show the resulting values for each metric. By default, the results are compared to the baseline. To generate this baseline, the cluster column of the clustering output is randomly permuted as to mimic a clustering algorithm that performs random clustering.

```python
print(retention)
>> 0.2853333333333333
print(purity)
>> {'Actual': 0.8916464891041163, 'Permuted': 0.21337772397094432}
print(consistency)
>> {'Actual': 0.4064769975786925, 'Permuted': 0.06174334140435835}
```
