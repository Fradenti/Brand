
<!-- README.md is generated from README.Rmd. Please edit that file -->

# brand

<!-- badges: start -->
<!-- badges: end -->

R package for fitting Brand: Bayesian Robust Adaptive Novelty Detector.
Brand is a two-stage Bayesian semiparametric model. In the first stage,
it learns the main characteristics of the known classes from the labeled
dataset using robust procedures. In the second phase, a Bayesian
semiparametric mixture of known groups and a novelty term is fitted to
the test set. Training insights is used to elicit informative priors for
the known components. The novelty term is instead captured via a
flexible Dirichlet Process mixture.

The package provides efficient Slice Samplers for handling both
multivariate and functional data.

The repository is associated with the paper Denti, Cappozzo, Greselin
(2021) *A two-stage Bayesian semiparametric model for novelty detection
with robust prior information.*
<https://doi.org/10.1007/s11222-021-10017-7>

For replicating the simulated and real data analyses reported in the
paper, please referer to
[this](https://github.com/AndreaCappozzo/brand-public_repo) repository.

## Installation

You can install the development version of brand from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("Fradenti/Brand")
```
