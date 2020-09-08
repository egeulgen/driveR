
<!-- README.md is generated from README.Rmd. Please edit that file -->

# <img src="https://github.com/egeulgen/driveR/blob/master/inst/extdata/driveR_logo.png?raw=true" align="left" height=150/> driveR: An R Package for Prioritizing Cancer Driver Genes Using Genomics Data

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/egeulgen/driveR.svg?branch=master)](https://travis-ci.com/egeulgen/driveR)
[![Codecov test
coverage](https://codecov.io/gh/egeulgen/driveR/branch/master/graph/badge.svg)](https://codecov.io/gh/egeulgen/driveR?branch=master)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

Cancer genomes contain large numbers of somatic alterations but few
genes drive tumor development. Identifying molecular cancer driver genes
is critical for precision oncology. Most of current approaches either
identify driver genes based on mutational recurrence or using estimated
scores predicting the functional consequences of mutations.

`driveR` is a tool for scoring and prioritizing genes from individual or
batch somatic whole exome/genome sequencing data according to cancer
driverness. As features, driveR uses coding impact metaprediction
scores, non-coding impact scores, somatic copy number alteration scores,
hotspot gene/double-hit gene condition, ‘phenolyzer’ gene scores and
memberships to cancer-related KEGG pathways. It uses these features to
estimate cancer-type-specific driverness probabilities for each gene
using the related task of a multi-task learning model.

## Installation

You can install the development version of `driveR` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("egeulgen/driveR", build_vignettes = TRUE)
```

## Usage

For detailed information on how to use `driveR`, please see the vignette
“How to use driveR” via `vignette("how_to_use")`
