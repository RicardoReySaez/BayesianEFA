# BayesianEFA (development version)

The v0.2.0 rewrite is in progress on the `dev` branch.

* The Stan models are no longer precompiled into the package. The package now
  installs in seconds and depends on Rcpp and RcppArmadillo only for its C++
  code. Stan sampling is unavailable until the new backends arrive: the models
  will be generated in R and compiled once per model.
* rstan, cmdstanr, stanr and stanli are optional backends (Suggests).
* `rsp_align()` rotates the draws with a C++ Varimax that repeats
  `stats::varimax(normalize = FALSE)` exactly, 4 to 12 times faster. The rest
  of the algorithm is unchanged.
* The pointwise log-likelihood of the null model and the EM estimate of the
  saturated moments with missing data run in C++.
* The Stan code is generated in R from the structure of the model (data
  model, prior on the loadings, patterns of missing values), with every
  additive constant, the normalizing constant of the unit-vector repulsion
  included. The likelihood needs one Cholesky decomposition per pattern of
  observed items and reads a root of the scatter matrix of each pattern, so
  FIML no longer evaluates every incomplete row on its own.
* Four backends (cmdstanr, stanli, stanr and rstan) run three methods
  (MCMC, ADVI and Pathfinder; rstan has no Pathfinder), chosen with
  `backend = "auto"` by default. Each generated model compiles once and is
  kept in `tools::R_user_dir("BayesianEFA", "cache")`; `befa_compile_models()`
  compiles them all in advance and `befa_clear_cache()` removes them. stanli
  needs no compilation.
* Every chain runs in its own background process, with its own starting
  values and seed, behind one progress bar that looks the same with every
  backend. No backend is loaded in the user's session.
* The chains start from the maximum-likelihood EFA that `stats::factanal()`
  and lavaan start from, with the unit-vector norms at the mode of the
  repulsion.
* The normal prior on the loadings now models `Sigma = Lambda Lambda' + Psi`
  directly. The location and scale of each loading follow the SD of its item,
  and the scale of the inverse-gamma prior on each uniqueness its variance, so
  the prior does not depend on the units of the items.

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
