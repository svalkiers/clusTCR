---
layout: default
title: Importing
nav_order: 3
---

##  Importing data

Data import functions are provided within the `datasets` method. A method for reading common rep-seq formats is also provided within clusTCR. Additionally, *clusTCR* offers a function that randomly samples sequences from a data set to construct metarepertoires.

```python
from clustcr import datasets, read_cdr3, metarepertoire
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

### Supported input formats

You can import rep-seq data from different sources using the `read_cdr3()` method.  *clusTCR* supports the following input formats.

| Format         | Method                                   | Info                                                         |
| -------------- | ---------------------------------------- | ------------------------------------------------------------ |
| immuneACCESS   | `read_cdr3(file, data_format='immuneaccess')` | Version automatically detected. More info about the immuneACCESS format: https://clients.adaptivebiotech.com/immuneaccess. |
| AIRR standards | `read_cdr3(file, data_format='airr')`         | More info about the AIRR standards data representation: https://docs.airr-community.org/en/stable/datarep/rearrangements.html. |
| TCRex          | `read_cdr3(file, data_format='tcrex')`        | More info about the TCRex format: https://tcrex.biodatamining.be/instructions/. |

Rep-seq files can be easily imported from the `.read_cdr3()` method. This data is not provided within the package, hence you will have to refer to the directory in which you have saved your files. Specify the input format and *clusTCR* will return a `pd.Series` of all unique CDR3 sequences from that file.

#### Importing a single rep-seq file

To parse and extract CDR3 sequences from a single file, you can use the `read_cdr3()` function.  This function takes two argument: the path to the file and the input format. For example, importing an immuneACCESS file might look something like this.

```python
data = read_cdr3('immuneACCESS_file.csv', data_format='immuneaccess')
```

#### Creating metarepertoires

*clusTCR* offers a functionality that can be used to create (large) metarepertoire files from a directory containing TCR repertoires. To create a metarepertoire, you must specificy the directory in which the files are located, as well as the input format, output format and desired size of the metarepertoire. *clusTCR* currently offers three output formats:

| Format  | Info                                                         |
| ------- | ------------------------------------------------------------ |
| CDR3    | Provides a `pandas.Series` of unique CDR3 sequences.         |
| GLIPH2  | Provides a `pandas.DataFrame` in the required input format for *GLIPH2*. |
| tcrdist | Provides a `pandas.DataFrame` in the required input format for *tcrdist3*. |

Here is an example of how this function can be used to sample one million unique CDR3 sequences from an immuneACCESS data set:

```python
meta = metarepertoire('immuneACCESS_file_directory/', data_format='immuneaccess', out_format='CDR3', n_sequences=10**6)
```

This will return a `pandas.Series` of (1 million) unique CDR3 sequences sampled from repertoires in the specified directory.
