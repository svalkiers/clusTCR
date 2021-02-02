---
layout: default
title: Importing data
nav_order: 3
---

### Built-in data

clusTCR provides a benchmark data set, which contains all human TRB entries from the VDJdb. Example data sets (which correspond to high-quality VDJdb entries) are available that can be used for benchmarking and or exploration.

```python
from clustcr.input.datasets import test_cdr3, test_epitopes
cdr3 = test_cdr3() # pd.Series of CDR3 sequences
epitope_data = test_epitopes() # CDR3 data with their corresponding epitopes
```

You can also import (a recent version of) the complete VDJdb:

```python
from clustcr.input.datasets import vdjdb_cdr3, vdjdb_epitopes
vdjdb_cdr3 = vdjdb_cdr3()
vdjdb_epitopes = vdjdb_epitopes()
```

### Importing immuneACCESS files

immuneACCESS files can be easily imported from the *datasets* module. This data is not provided within the package, hence you will have to refer to the directory in which you have saved your immuneACCESS files. clusTCR automatically detects the version of the immuneACCESS file (v1 or v2) and parses the file accordingly.

#### Importing a single repertoire file in the immuneACCESS format

To parse and extract CDR3 sequences from a single immuneACCESS file, you can use the `immuneACCESS_cdr3()` function.

```python
from clustcr.input.datasets import immuneACCESS_cdr3
immuneaccess_repertoire = immuneACCESS_cdr3()
```

#### Creating metarepertoires

clusTCR offers a functionality that can be used to create (large) metarepertoire files from a directory containing TCR repertoires in the immuneACCESS format. To create a metarepertoire, you must specificy the directory in which the immuneACCESS files are located, as well as the desired size of the metarepertoire.

```python
from clustcr.input.datasets import metarepertoire_cdr3
metarepertoire = metarepertoire_cdr3('immuneACCESS_file_directory/', n_sequences=10**6)
```

This will return a list of unique CDR3 sequences sampled from repertoires in the specified directory.

