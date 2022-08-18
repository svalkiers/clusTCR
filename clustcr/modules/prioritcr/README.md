# A Python toolkit for interactively exploring and dissecting T-cell receptor repertoires.

<img src="examples/example.png" alt="J" width="80%"/>

## Introduction

High-throughput T-cell receptor sequencing allows for the comprehensive
characterization of T-cell repertoires; yet, extracting their useful or
clinically actionable information remains particularly challenging, largely
owing to their extreme complexity and diversity. Recent methods that cluster
TCR CDR3 sequences provide a promising approach to reduce repertoire complexity
while retaining this information. This package aims to provide a framework with
diverse tools and functions for advanced and efficient post-analysis of such
clustered repertoires. Overall, this functionality comprises two components:

1. Different functions that define similarity between clusters and locality
sensitive hashing implementations that allow quick retrieval of similar
clusters. Together, these functions and methods form the basis of the
`ClusterRepertoireVisualization` class, providing an interactive, interpretable
yet holistic overview of the cluster repertoire.

2. An extensive range of functions to efficiently calculate cluster features,
ranging from e.g. their generation probability to cluster convergence.
Together, these features allow a deep dive into the TCR immune response, and
provide means to filter clusters by features that might indicate antigen-driven
convergent recombination – possibly providing researchers with a method to
prioritize highly informative TCRs from complex repertoires. 

## Examples
