# TCR Clustering

This project aims to efficiently and accurately cluster T-cell receptors using different distance metrics.  
It mainly relies on the [faiss](https://github.com/facebookresearch/faiss/) library to cluster the sequences.

This project is made during an Internship at the [Adrem Data Lab](https://www.uantwerpen.be/en/research-groups/adrem-data-lab/) 
of the University of Antwerp with supervision and help of [Pieter Meysman](https://www.uantwerpen.be/nl/personeel/pieter-meysman/) and others.

 
 ## Install
 
 The first step to be able to run this project is installing the faiss library.
 For this they have an [instruction manual](https://github.com/facebookresearch/faiss/blob/master/INSTALL.md).
 The easiest way to install is using conda which comes down to
 
 ```
 conda install faiss-cpu -c pytorch
 ```
 
 If the installation fails because of an incorrect python version, make a new conda python environment with a version that is supported.
 This should fix the problem.
 
 
 
 
