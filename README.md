
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
#> x5          0.86*                0.75  0.25  1.00     3716     3149
#> x4          0.83*                0.71  0.29  1.00     3756     3125
#> x6          0.81*                0.69  0.31  1.00     3943     3132
#> x8                 0.69*         0.53  0.47  1.00     2655     2169
#> x7                 0.68*         0.49  0.51  1.00     2864     2196
#> x9                 0.49*  0.40*  0.44  0.56  1.00     3834     2732
#> x3                        0.65*  0.45  0.55  1.00     3810     2267
#> x1          0.31*         0.60*  0.48  0.52  1.00     3744     2880
#> x2                        0.47*  0.25  0.75  1.00     3212     2673
#> ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗ 
#> Note: varimax rotation applied; MCMC with stanli. Diagnostics show
#> worst-case values across factors (max Rhat, min ESS). (*) 95%
#> Credible Interval excludes 0. Loadings with absolute values < 0.30
#> are hidden. 
#> 
#> Table 2. Variance Explained 
#> ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗ 
#> Statistic                    F1    F2    F3  Total
#> ——————————————————————————————————————————————————
#> Sum of squared loadings    2.27  1.27  1.26   4.79
#> Proportion of total        0.47  0.26  0.26   1.00
#> Proportion of variance     0.25  0.14  0.14   0.53
#> Cumulative variance        0.25  0.39  0.53       
#> ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗ 
#> Note: Posterior means. Proportion of total: share
#> of the explained variance. 
#> 
#> Table 3. Bayesian Fit Measures 
#> ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗ 
#> Index         Estimate    SD  CI_Low  CI_High
#> —————————————————————————————————————————————
#> Chi2             47.05  6.92   35.44    62.49
#> Chi2_ppp          0.11                       
#> Chi2_Null       918.85  0.00  918.85   918.85
#> BRMSEA            0.05  0.02    0.00     0.09
#> BGamma            0.99  0.00    0.98     1.00
#> Adj_BGamma        0.97  0.02    0.94     1.00
#> BMc               0.98  0.01    0.96     1.00
#> SRMR              0.05  0.01    0.03     0.06
#> BCFI              0.99  0.01    0.97     1.00
#> BTLI              0.99  0.02    0.94     1.00
#> ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗ 
#> Note: Intervals are 95% credible intervals.
#> PPP: Posterior Predictive p-value (Ideal >
#> .05). 
#> 
#> Table 4. Information Criteria 
#> ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗ 
#> Criterion    Estimate     SE    CI_Low   CI_High
#> ————————————————————————————————————————————————
#> ELPD         -3420.92  42.68  -3504.58  -3337.27
#> LOOIC         6841.84  85.36   6674.54   7009.15
#> p_loo           25.26   1.79     21.76     28.76
#> WAIC          6841.75  85.35   6674.47   7009.04
#> p_WAIC          25.21   1.78     21.72     28.70
#> DIC           6839.65  84.04   6674.93   7004.37
#> p_D             23.64   0.80     22.08     25.21
#> DIC_i         6839.94  79.98   6683.17   6996.70
#> p_V             23.93   5.57     13.03     34.84
#> BPIC          6863.29  84.41   6697.85   7028.73
#> ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗ 
#> Note: SEs across persons (jackknife for DIC_i
#> and p_V); 95% normal intervals. 
#> 
#> Table 5. Factor Reliability (Omega) 
#> ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗ 
#> Factor    Estimate    SD  CI_Low  CI_High
#> —————————————————————————————————————————
#> F1            0.73  0.02    0.69     0.76
#> F2            0.55  0.04    0.45     0.62
#> F3            0.60  0.04    0.52     0.67
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
