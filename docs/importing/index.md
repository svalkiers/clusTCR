---
layout: default
title: Importing
nav_order: 3
---

##  Importing data

Data import functions are provided within the `datasets` method.

```python
from clustcr import datasets
```

### Built-in data

clusTCR provides a benchmark data set, which contains all human TRB entries from the VDJdb. Example data sets (which correspond to high-quality VDJdb entries) are available that can be used for benchmarking and or exploration.

```python
cdr3 = datasets.test_cdr3() # pd.Series of CDR3 sequences
epitope_data = datasets.test_epitopes() # CDR3 data with their corresponding epitopes
```

You can also import (a recent version of) the complete VDJdb:

```python
vdjdb_cdr3 = datasets.vdjdb_cdr3()
vdjdb_epitopes = datasets.vdjdb_epitopes()
```

### Importing immuneACCESS files

immuneACCESS files can be easily imported from the `datasets` module. This data is not provided within the package, hence you will have to refer to the directory in which you have saved your immuneACCESS files. clusTCR automatically detects the version of the immuneACCESS file (v1 or v2) and parses the file accordingly.

#### Importing a single repertoire file in the immuneACCESS format

To parse and extract CDR3 sequences from a single immuneACCESS file, you can use the `immuneACCESS_cdr3()` function.  This function takes one argument: the path to the immuneACCESS file.

```python
immuneaccess_repertoire = datasets.immuneACCESS_cdr3('immuneACCESS_file.csv')
```

#### Creating metarepertoires

clusTCR offers a functionality that can be used to create (large) metarepertoire files from a directory containing TCR repertoires in the immuneACCESS format. To create a metarepertoire, you must specificy the directory in which the immuneACCESS files are located, as well as the desired size of the metarepertoire.

```python
metarepertoire = datasets.metarepertoire_cdr3('immuneACCESS_file_directory/', n_sequences=10**6)
```

This will return a list of unique CDR3 sequences sampled from repertoires in the specified directory.

