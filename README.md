# mavric

MAVRIC quantifies the percentage of variance due to each experimental covariate in genomics data. The method further provides estimates of confounding effects and nonparametric statistical testing of differential counts. Details are provided in a preprint available on [bioRxiv](https://www.biorxiv.org/content/early/2018/05/04/314112), with an associated implementation in R provided here.

## Usage

The primary user-facing function in the package is `estVC`, which performs all the calculations underlying the algorithm and returns the above-described quantification and differential analysis. Additionally, `plotPCs` and `plotVars` offer visualizations of the results. All three functions have help pages within R describing their parameters and use in greater detail.
