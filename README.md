# mavric

MAVRIC quantifies the percentage of variance due to each experimental covariate in genomics data. The method further provides estimates of confounding effects and nonparametric statistical testing of differential counts. Details are provided in a preprint available on [bioRxiv](https://www.biorxiv.org/content/early/2018/05/04/314112), with an associated implementation in R provided here.

## Installation

MAVRIC can be loaded into R using the `load_all` function in the [devtools](https://cran.r-project.org/web/packages/devtools/index.html) R package, directing the function at the base directory of the package. The requirements for MAVRIC are enumerated in the DESCRIPTION file; once those dependencies have been installed within R, the package should load within seconds. MAVRIC has been evaluated on Ubuntu 14.04 (Trusty Tahr).

## Usage

The primary user-facing function in the package is `estVC`, which performs all the calculations underlying the algorithm and returns the above-described quantification and differential analysis. Additionally, `plotPCs` and `plotVars` offer visualizations of the results. All three functions have help pages within R describing their parameters and use in greater detail. The help page for `estVC` includes an example for running all three functions.
