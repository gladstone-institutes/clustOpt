# clustOpt <img src="man/figures/clustOpt_logo.png" align="right" height="138" alt="clustOpt logo" /></a>

<!-- badges: start -->
[![R-CMD-check](https://github.com/gladstone-institutes/clustOpt/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gladstone-institutes/clustOpt/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->
   

*Authors: Natalie Gill, Min-Gyoung Shin,Ayushi Agrawal, and Reuben Thomas*


Selecting the clustering resolution parameter for Louvain clustering in scRNA-seq is often based on the concentration of expression of cell type marker genes within clusters, increasing the parameter as needed to resolve clusters with mixed cell type gene signatures. This approach is however subjective in situations where one does not have complete knowledge of condition/disease associated cell-types in the context of novel biology, it is time-consuming and has the potential to bias the final clustering results due to individual transcriptomic heterogeneity, and subject-specific differences in cell composition.

clustOpt improves the reproducibility of modularity based clustering in multi-subject experiments by using a combination of subject-wise cross validation, feature splitting, random forests and measures of cluster quality using the silhouette metric to guide the selection of the resolution parameter. 

## clustOpt Algorithm

 <img src="man/figures/clustOpt_diagram.png" align="center" alt="clustOpt algorithm" /></a>


To avoid the issue of data leakage we cluster the evaluate in the odd PC space while training and predicting on the holdout subject in the even PC space.

<p align="center">
<img src="man/figures/pc_split.png" width="50%" alt="clustOpt logo" />
</p>

## Installation   
Currently the only way to install is by using the package`devtools`:    
```
devtools::install_github("gladstone-institutes/clustOpt")
```
If you get an error message and everything is spelled correctly, follow these steps before trying again:
```
#set config
usethis::use_git_config(user.name = "YourName", user.email = "your@mail.com")

#Go to github page to generate token
usethis::create_github_token() 

#paste your PAT into pop-up that follows...
credentials::set_github_pat()
```

## Get Started

