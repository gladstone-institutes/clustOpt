<!-- badges: start -->
[![R-CMD-check](https://github.com/gladstone-institutes/clustOpt/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gladstone-institutes/clustOpt/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->
   
# clustOpt
*Authors: Natalie Elphick, Min-Gyoung Shin,Ayushi Agrawal, and Reuben Thomas*

This package provides a way to choose a resolution parameter that minimizes subject-specific effects on cluster membership. It uses random forest (RF) models with an approach similar to leave one out cross validation, where each subject is held out and the rest are used to train a RF to predict the cluster labels for the held out subject. This is done across a range of resolution values to generate a distribution of silhouette scores for the held out subject across the iterations. This silhouette score distribution is then used to choose a resolution value whose mean silhouette score results in a local maximum.

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

# Get Started

