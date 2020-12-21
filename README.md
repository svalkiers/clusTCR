# TCR clustering pipeline
Two-step TCR clustering approach that combines the speed of the [faiss](https://github.com/facebookresearch/faiss) library with the accuracy of the [Markov Clustering Algorithm](https://micans.org/mcl/).

## Dependencies
The minimum requirement to use `networkTCR.py` is the installation of the following Python packages:
<br/>
- Through `conda`
```
conda install pandas
conda install networkx
conda install scikit-learn
conda install faiss-cpu -c pytorch
```
- Through `pip`
```
pip install markov_clustering
pip install pyteomics
```