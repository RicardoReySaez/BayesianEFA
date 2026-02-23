
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BayesEFA: Bayesian Exploratory Factor Analysis

<!-- badges: start -->

<!-- badges: end -->

The goal of **BayesEFA** is to provide a flexible and robust framework
for conducting Bayesian Exploratory Factor Analysis (EFA). While
remaining as intuitive and easy to use as standard frequentist software,
**BayesEFA** offers the significant advantage of solving common
estimation issues by default, such as non-positive definite correlation
matrices and Heywood cases.

## Installation

### 1. Install `cmdstanr` (Recommended)

While **BayesEFA** can use the `rstan` backend, **we strongly recommend
using the `cmdstanr` backend** for faster estimation and lower memory
usage.

Please install `cmdstanr` by following the official guide:

> [**Getting Started with
> CmdStanR**](https://mc-stan.org/cmdstanr/articles/cmdstanr.html)

### 2. Install BayesEFA

You can install the development version of BayesEFA from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("https://github.com/RicardoReySaez/BayesEFA")
```

Once installed, **BayesEFA** will automatically detect and use
`cmdstanr` if available.

## Getting Started

For a complete guide on how to use the package, please visit the [Get
Started](articles/getting_started.html) page or the package website.
