---
layout: default
title: Importing
nav_order: 3
---

##  Importing data

Data import functions are provided within the `datasets` method. A method for reading common rep-seq formats is also provided within clusTCR. Additionally, ClusTCR offers a function that randomly samples sequences from a data set to construct metarepertoires.

```python
from clustcr import datasets, read_cdr3, metarepertoire
```

### Built-in data

#### Test data

ClusTCR provides a test data set, which contains all high-quality human TCRB entries from the VDJdb. The example data sets can be used for benchmarking and/or exploration.

```python
cdr3 = datasets.test_cdr3() # pd.Series of CDR3 sequences
epitope_data = datasets.test_epitopes() # CDR3 data with their corresponding epitopes
```

#### VDJdb

A complete, recent version of the VDJdb is provided within the package. Users can specify whether they want to solely import beta, alpha or paired chains. Corresponding epitope information can also be retrieved.

##### TCR sequences

You can import alpha, beta or paired chain sequences from VDJdb:

```python
vdjdb_b = datasets.vdjdb_beta() # beta chain
vdjdb_a = datasets.vdjdb_alpha() # alpha chain
vdjdb_ab = datasets.vdjdb_paired() # paired alpha-beta

print(vdjdb_b.head()) # example
```

This will return a `pandas.Series` of TCRB sequences:

```
0    CASSYLPGQGDHYSNQPQHF
1     CASSFEAGQGFFSNQPQHF
2     CASSFEPGQGFYSNQPQHF
3    CASSYEPGQVSHYSNQPQHF
4           CASSFGVEDEQYF
Name: CDR3_beta, dtype: object
```

##### Epitope information

Similarly, you can add epitope information by specifying the `epitopes=True` argument.

```python
vdjdb_b_epitopes = datasets.vdjdb_beta(epitopes=True)
print(vdjdb_b_epitopes.head())
```

This will return a two-column `pandas.DataFrame`:

```
              CDR3_beta   Epitope
0  CASSYLPGQGDHYSNQPQHF  FLKEKGGL
1   CASSFEAGQGFFSNQPQHF  FLKEKGGL
2   CASSFEPGQGFYSNQPQHF  FLKEKGGL
3  CASSYEPGQVSHYSNQPQHF  FLKEKGGL
4  CASSYLPGQGDHYSNQPQHF  FLKEQGGL
```

##### Quality cut-off

All VDJdb entries are assigned a score (*q*) that represents the confidence of the recorded TCR-epitope interaction. The VDJdb data incorporated in ClusTCR can be filtered based on the *q* score by specifying a value for the `q` argument. The authors of VDJdb give the following descriptions for the different scores:

| **score** | description                                                  |
| --------- | ------------------------------------------------------------ |
| 0         | Low confidence/no information - a critical aspect of sequencing/specificity validation is missing |
| 1         | Moderate confidence - no verification / poor TCR sequence confidence |
| 2         | High confidence - has some specificity verification, good TCR sequence confidence |
| 3         | Very high confidence - has extensive verification or structural data |

Example:

```python
vdjdb_b_highqual = datasets.vdjdb_beta(q=2)
```

### Importing a repertoire file

ClusTCR groups TCR sequences **based on their CDR3 amino acid sequence**, V/J gene information is not used. To parse and extract CDR3 sequences from a repertoire sequencing (rep-seq) file, you can use the `.read_cdr3()` method. This function takes two argument: the path to the file and the input format. For example, importing an immuneACCESS file might look something like this:

```python
data = read_cdr3('immuneACCESS_file.csv', data_format='immuneaccess')
```

ClusTCR will return a `pandas.Series` of all unique CDR3 sequences from that file. 

#### Supported input formats

You can import rep-seq data from different sources using the `read_cdr3()` method.  *clusTCR* supports the following input formats.

| Format         | Method                                        | Info                                                         |
| -------------- | --------------------------------------------- | ------------------------------------------------------------ |
| immuneACCESS   | `read_cdr3(file, data_format='immuneaccess')` | Version automatically detected. More info about the immuneACCESS format: https://clients.adaptivebiotech.com/immuneaccess. |
| AIRR standards | `read_cdr3(file, data_format='airr')`         | More info about the AIRR standards data representation: https://docs.airr-community.org/en/stable/datarep/rearrangements.html. |
| TCRex          | `read_cdr3(file, data_format='tcrex')`        | More info about the TCRex format: https://tcrex.biodatamining.be/instructions/. |

