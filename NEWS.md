# BayesianEFA (development version)

The v0.2.0 rewrite is in progress on the `dev` branch.

* The Stan models are no longer precompiled into the package. The package now
  installs in seconds and depends on Rcpp and RcppArmadillo only for its C++
  code. Stan sampling is unavailable until the new backends arrive: the models
  will be generated in R and compiled once per model.
* rstan, cmdstanr, stanr and stanli are optional backends (Suggests).

# BayesianEFA 0.1.0

* Initial CRAN release.

## Features

* Bayesian Exploratory Factor Analysis via Hamiltonian Monte Carlo (Stan).
* Three model parameterizations: correlation (`"cor"`), covariance (`"cov"`),
  and full unstandardized (`"raw"`).
* Unit-vector prior for factor loadings ensuring valid variance partitioning
  and Heywood case prevention (Rey-Sáez et al., 2025).
* RSP alignment algorithm with Varimax rotation for resolving rotational
  indeterminacy across MCMC draws.
* Full Information Maximum Likelihood (FIML) for missing data handling.
* Bayesian fit indices: BRMSEA, BCFI, BTLI, SRMR, Chi-squared with
  posterior predictive p-value, and PSIS-LOO cross-validation.
* Omega reliability coefficients (total and subscale) with full posterior
  uncertainty quantification.
* Bayesian factor scores with posterior distributions.
* Dual backend support: `rstan` and `cmdstanr`.
* Efficient C++ implementations (RcppArmadillo) for the RSP algorithm
  and factor score computation.
