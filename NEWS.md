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
* `befa()` runs on the new pipeline. New arguments `method` (`"mcmc"`,
  `"advi"`, `"pathfinder"`) and `backend` (default `"auto"`); MCMC options use
  the cmdstanr names (`chains`, `parallel_chains`, `iter_warmup`,
  `iter_sampling`, `adapt_delta`, `max_treedepth`, ...) and still accept the
  rstan ones. Unknown options are errors.
* Breaking: the fit no longer holds a `stanfit`. `$draws` holds the aligned
  parameters in the metric of the model (`Lambda_std`, `Psi_std`, `h2` for
  the correlation model; `Lambda`, `Psi`, `sigma`, `h2`, `nu` otherwise), and
  `extract_posterior_draws()` and `posterior_summaries()` gain
  `metric = c("model", "std", "raw")`.
* After the fit, `befa()` reports divergent transitions, maximum tree depth,
  the largest R-hat (below 1.01) and the smallest effective sample size
  (above 400), computed on the aligned parameters only, and warns with the
  remedy when a check fails. ADVI and Pathfinder report the Pareto k of the
  approximation.
* Factor scores, reliability, fit indices and `summary()` read the final
  (aligned) draws. `befa_factor_scores()`, `befa_reliability()`,
  `befa_fit_measures()` and `summary()` work on any fit, and `befa()` stores
  their results unless asked not to.
  * The factor scores are drawn in the metric the model was fitted to,
    follow `seed`, and are also returned by `extract_posterior_draws()`.
  * Omega uses the standardized solution in every data model.
  * The null model of the incremental fit indices (BCFI, BTLI) is fitted
    with the backend and method of the fit; the null correlation model has
    no free parameters and needs no fit. PSIS-LOO treats the draws of ADVI
    and Pathfinder as independent.
  * `summary()` shows the standardized loadings, communalities and
    uniquenesses in every data model, with the convergence of each item and
    the problems the fit reported.
  * Removed: the internal helpers that read a `stanfit`.
* The normal prior on the loadings now models `Sigma = Lambda Lambda' + Psi`
  directly. The location and scale of each loading follow the SD of its item,
  and the scale of the inverse-gamma prior on each uniqueness its variance, so
  the prior does not depend on the units of the items.
* The likelihood uses the Woodbury identity, which works with an M x M
  matrix instead of the J x J covariance matrix. One gradient costs up to 14
  times less with complete data, and under FIML up to 175 times less, as
  everything that grows with the number of items is computed for all the
  patterns of missing values at once. The benchmark is in `bench/`, outside
  the package.
* New defaults: `backend = "auto"` tries stanli first (then cmdstanr, stanr
  and rstan), since it needs no compilation and samples these models as fast
  as the compiled backends; MCMC uses `adapt_delta = 0.9` (was 0.95), which
  gave more effective draws per second without divergent transitions.
* Pathfinder with cmdstanr or stanr reports a finite Pareto k: the PSIS
  resampling now happens after k is computed on all the candidate draws.
* The vignettes and the README use the new interface. The vignettes are
  precomputed (`vignettes/*.Rmd.orig`, knitted by `vignettes/precompile.R`),
  so building the package and R CMD check run no Stan. The README shows the
  Codecov badge, and the pkgdown reference is grouped by topic.
* The aligned factors are ordered by the variance they explain, so F1 is the
  strongest factor in the draws, `summary()`, the factor scores and omega.
* `summary()` shows the names of the items (from the columns of `data` or of
  the summary statistics) instead of `Item_1`, ..., in a column of fixed
  width (`name_width`, default 8) that keeps the width of the table: longer
  names are shortened in the middle. The items are sorted by their main
  factor by default (`sort = TRUE`), omega total is the last row of the
  reliability table, and every note wraps to the width of its table.
* Information criteria: `befa_fit_measures()` adds WAIC, the DIC of
  Spiegelhalter et al. (2002) with the posterior means of the aligned
  parameters, the parameterization-invariant DIC_i of Xiao and Rabe-Hesketh
  (2026), and the BPIC, to ELPD, LOOIC and p_loo, each with its standard
  error across persons (by the jackknife for DIC_i and p_V). `summary()`
  shows them in their own table, between the fit measures and the
  reliability.
* New `befa_bridge_sampler()`: the log marginal likelihood by bridge
  sampling (through the bridgesampling package), for Bayes factors between
  numbers of factors. Since the posterior is invariant under rotations of
  the loadings, every draw is rotated to a canonical representative (QR
  decomposition of a block of anchor rows) and the integral over the
  orthogonal group is taken exactly, so the raw draws need not visit every
  equivalent copy. The log density runs in C++ with the vectorized Woodbury
  likelihood of the Stan code. `befa()` now keeps the draws of the Stan
  parameters in `$fit$draws`, which this needs.
* New `befa_compare()`: compares two or more fits to the same data in four
  tables, each ordered from the best model: the posterior means of the
  Bayesian fit indices (by BRMSEA), the PSIS-LOO comparison (differences
  from the best model, `p_worse`, and stacking and pseudo-BMA+ weights), the
  information criteria with Akaike weights (by WAIC), and the model evidence
  from bridge sampling with posterior model probabilities. The Monte Carlo
  error of each evidence is the closed form of Micaletto and Vehtari (2025)
  from bridgesampling 1.2.1, so no repetitions are needed; the probability
  intervals propagate it, and prior model probabilities can be given.
* `befa()` accepts several values of `n_factors`, such as `n_factors = 1:3`:
  it fits one model per value with the same settings and seed, and returns a
  `befa_list` with the fits as `nf1`, `nf2`, ... The saturated and null models
  are estimated once. `compare = TRUE` (the default) runs `befa_compare()` with
  `compare_criteria = c("fit", "loo", "ic")`; add `"evidence"` for the posterior
  model probabilities. A model that fails does not stop the others.
  `summary()` shows an overview of the models, the comparison tables and the
  loadings of the model with the highest ELPD (`show` chooses others).
  `befa_compare()`, `befa_fit_measures()` and `befa_reliability()` accept the
  list.
* `summary()` of a fit adds a table of variance explained (sum of squared
  loadings and proportions), computed over the draws, so its total equals the
  sum of the communalities. Factors are ordered by the same quantity.
* The references cite the published papers: the standardized generalized
  hierarchical factor model in *Psychological Methods*
  (doi:10.1037/met0000846) and the E-RSP algorithm in *Bayesian Analysis*
  (doi:10.1214/26-BA1611).

# BayesianEFA 0.1.0

* Initial CRAN release.

## Features

* Bayesian Exploratory Factor Analysis via Hamiltonian Monte Carlo (Stan).
* Three model parameterizations: correlation (`"cor"`), covariance (`"cov"`),
  and full unstandardized (`"raw"`).
* Unit-vector prior for factor loadings ensuring valid variance partitioning
  and Heywood case prevention (Rey-Sáez et al., 2026).
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
