
<!-- README.md is generated from README.Rmd. Please edit that file -->

An R package for Effect Consistency Test

## Installation

You can install the development version of ect from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("xsun1229/ect")
```

Then load it

``` r
library(ect)
```

## Overview of the Functions

We provide two main functions in this package:

1. `prune_snps()` — Select independent SNPs

Selects independent SNPs using PLINK clumping and LD pruning.

2. `ECT()` — Effect Consistency Test

Performs an adaptive resampling test to assess whether SNP effect directions on an outcome (GWAS trait) are consistent with those on an exposure (pathway factor).


## Full tutorial

<https://xsun1229.github.io/ect/articles/getting_started.html>

<!-- ## Citation -->
<!-- If you use ect in published work, please cite it using: -->
<!-- ``` r -->
<!-- citation("ect") -->
<!-- ``` -->
