
<!-- README.md is generated from README.Rmd. Please edit that file and run
     vignettes/precompile.R, which also knits the vignettes. -->

# BayesianEFA: Bayesian Exploratory Factor Analysis <img src="man/figures/logo.png" align="right" height="120" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/RicardoReySaez/BayesianEFA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/RicardoReySaez/BayesianEFA/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/RicardoReySaez/BayesianEFA/graph/badge.svg)](https://app.codecov.io/gh/RicardoReySaez/BayesianEFA)
[![CRAN
status](https://www.r-pkg.org/badges/version/BayesianEFA)](https://CRAN.R-project.org/package=BayesianEFA)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

The `BayesianEFA` package estimates Bayesian Exploratory Factor Analysis
(EFA) models with `Stan`. The loading matrix is left unrestricted, and
its rotational indeterminacy is solved after sampling with the Efficient
Rotation-Sign-Permutation (E-RSP) alignment. The package also provides
Bayesian fit indices, omega reliability coefficients and factor scores,
each with its full posterior distribution.

## Installation

You can install the development version of BayesianEFA from
[GitHub](https://github.com/RicardoReySaez/BayesianEFA) with:

``` r
# install.packages("devtools")
devtools::install_github("RicardoReySaez/BayesianEFA")
```

The package installs in seconds: its Stan programs are generated when
they are needed, not compiled during installation. Its own C++ code
needs the usual R toolchain (Rtools on Windows; see below for macOS).

## Stan backends

`BayesianEFA` runs Stan through one of four optional interfaces. We
recommend **stanli**, which interprets Stan programs, needs no C++
compiler and starts at once:

``` r
install.packages("stanli", repos = c("https://seantalts.r-universe.dev",
                                     getOption("repos")))
stanli::stanli_install()
```

[cmdstanr](https://mc-stan.org/cmdstanr/), stanr and rstan work too.
They compile each model the first time it is needed and keep it in a
cache for later sessions; `befa_compile_models()` compiles them all in
advance. With `backend = "auto"` (the default), `befa()` uses the first
one installed among stanli, cmdstanr, stanr and rstan.

## Installation in macOS

Since `BayesianEFA` includes custom C++ functions, macOS users must
configure the R compilation toolchain. Following the recommendations
from the [latent](https://github.com/Marcosjnez/latent) package, we
suggest using `macrtools` by James Balamuta to automate this setup:

``` r
# 1. Install macrtools
# install.packages("remotes")
remotes::install_github("coatless-mac/macrtools")

# 2. Setup the toolchain and OpenMP support
macrtools::macos_rtools_install()
macrtools::openmp_install()

# 3. Install BayesianEFA
devtools::install_github("RicardoReySaez/BayesianEFA")
```

## Quick Start

Here is a basic example showing how to fit a Bayesian EFA model:

``` r
library(BayesianEFA)

# Fit a 3-factor Bayesian EFA model to the Holzinger & Swineford data
befa_fit <- befa(data = HS_data, n_factors = 3, seed = 17)

# Print a summary of the model (inspired by the output of the psych package)
summary(befa_fit, cutoff = 0.3, signif_stars = TRUE)
#> 
#> Table 1. Standardized Factor Loadings 
#> ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗ 
#> Variable       F1     F2     F3    h2    u2  Rhat  EssBulk  EssTail
#> ———————————————————————————————————————————————————————————————————
#> x5          0.86*                0.75  0.25  1.00     3919     3085
#> x4          0.83*                0.72  0.28  1.00     3815     2080
#> x6          0.81*                0.69  0.31  1.00     4204     3187
#> x8                 0.70*         0.53  0.47  1.00     3074     1880
#> x7                 0.68*         0.48  0.52  1.00     3400     2490
#> x9                 0.49*  0.40*  0.44  0.56  1.00     3271     3014
#> x3                        0.65*  0.45  0.55  1.00     4054     3045
#> x1          0.31*         0.59*  0.48  0.52  1.00     3888     3205
#> x2                        0.46*  0.24  0.76  1.00     3846     2339
#> ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗ 
#> Note: varimax rotation applied; MCMC with stanli. Diagnostics show
#> worst-case values across factors (max Rhat, min ESS). The 3 latent
#> factors accounted for 52.2% of total variance. (*) 95% Credible
#> Interval excludes 0. Loadings with absolute values < 0.30 are
#> hidden. 
#> 
#> Table 2. Bayesian Fit Measures 
#> ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗ 
#> Index         Estimate    SD  CI_Low  CI_High
#> —————————————————————————————————————————————
#> Chi2             47.50  6.98   35.80    63.10
#> Chi2_ppp          0.11                       
#> Chi2_Null       918.85  0.00  918.85   918.85
#> BRMSEA            0.06  0.02    0.00     0.09
#> BGamma            0.99  0.00    0.98     1.00
#> Adj_BGamma        0.97  0.02    0.93     1.00
#> BMc               0.98  0.01    0.96     1.00
#> SRMR              0.05  0.01    0.03     0.06
#> BCFI              0.99  0.01    0.97     1.00
#> BTLI              0.99  0.02    0.93     1.00
#> ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗ 
#> Note: Intervals are 95% credible intervals.
#> PPP: Posterior Predictive p-value (Ideal >
#> .05). 
#> 
#> Table 3. Information Criteria 
#> ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗ 
#> Criterion    Estimate     SE    CI_Low   CI_High
#> ————————————————————————————————————————————————
#> ELPD         -3421.23  42.63  -3504.79  -3337.66
#> LOOIC         6842.45  85.27   6675.33   7009.57
#> p_loo           25.44   1.80     21.91     28.97
#> WAIC          6842.37  85.26   6675.27   7009.48
#> p_WAIC          25.40   1.80     21.88     28.92
#> DIC           6840.38  83.92   6675.90   7004.85
#> p_D             23.92   0.82     22.31     25.52
#> DIC_i         6840.83  79.61   6684.79   6996.87
#> p_V             24.37   6.03     12.56     36.18
#> BPIC          6864.29  84.26   6699.14   7029.45
#> ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗ 
#> Note: SEs across persons (jackknife for DIC_i
#> and p_V); 95% normal intervals. 
#> 
#> Table 4. Factor Reliability (Omega) 
#> ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗ 
#> Factor    Estimate    SD  CI_Low  CI_High
#> —————————————————————————————————————————
#> F1            0.73  0.02    0.69     0.76
#> F2            0.54  0.04    0.45     0.62
#> F3            0.60  0.04    0.51     0.67
#> Total         0.84  0.01    0.82     0.86
#> ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗ 
#> Note: Total is the omega of the full
#> scale. Omega coefficients use the
#> standardized solution.
```

For a complete guide on how to perform a Bayesian EFA from scratch,
please visit the [Getting
Started](https://ricardoreysaez.github.io/BayesianEFA/articles/getting_started.html)
tutorial or explore the package website.
