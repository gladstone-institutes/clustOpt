<!-- badges: start -->
[![R-CMD-check](https://github.com/gladstone-institutes/clustOpt/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gladstone-institutes/clustOpt/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->
   
# clustOpt
*Authors: Min-Gyoung Shin, Reuben Thomas, Natalie Elphick and Ayushi Agrawal*

## Background
Seurat uses a modularity optimization based clustering algorithm to cluster cells in a shared nearest neighbor graph. By default it uses a C++ implementation of the Louvian algorithm ([Blondel et al. 2008](https://doi.org/10.1088/1742-5468/2008/10/P10008)) described in [A smart local moving algorithm for large-scale modularity-based community detection](http://www.ludowaltman.nl/slm/) by Waltman and van Eck (2013). It first calculates the k-nearest neighbors, construct the SNN graph, and then optimizes the modularity function to determine clusters.

The resolution parameter controls the granularity of the communities (clusters) that are detected, larger resolution values result in more clusters. Choosing the correct resolution parameter for a given scRNA-seq dataset should result in clusters defined by distinct transcriptomic states such as cell types. Researchers typically use marker genes and the stability of the clusters that express them to choose the resolution parameter that best separates known cell types into their own clusters.

There are many potential sources of noise in scRNA-seq data that can result in clusters that are technical artifacts. The stability of a given cluster can also be impacted by a sample-specific variation. Since downstream analysis relies heavily on the assumption that clusters are defined by the biological variation of interest, it is important to choose a resolution that minimizes the effect of technical variation on cluster membership.

## Goal
This package provides a way to choose a resolution paramter that minimizes sample-specific effects on cluster membership. It uses random forest (RF) models with an approach similar to leave one out cross validation, where each sample is held out and the rest are used to train a RF to predict the cluster labels for the held out sample. This is done across a range of resolution values to generate a distribution of silhouette scores for the held out sample across the iterations. This silhouette score distribution is then used to choose a resultion value whose mean silhouette score results in a local maximum.

## Current Status
Currently implemented in an Rscript, this repo will be used to turn it into a package and incorporate some performance improvements.


