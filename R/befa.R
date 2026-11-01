# befa.R
# Fit a Bayesian exploratory factor analysis
# Author: Ricardo Rey-Sáez
# Last modified: 01-11-2026

# Public functions

#' Bayesian Exploratory Factor Analysis
#'
#' Fits an exploratory factor model with Stan and returns posterior draws of
#' interpretable parameters. The loading matrix is left unrestricted during
#' estimation; its rotational indeterminacy is resolved afterwards, draw by
#' draw, with Varimax and the Efficient Rotation-Sign-Permutation (E-RSP)
#' alignment of [rsp_align()].
#'
#' @param data Numeric matrix or data frame with one row per person and one
#'   column per item. Leave it `NULL` to fit summary statistics instead.
#' @param n_factors Number of latent factors, smaller than the number of
#'   items.
#' @param ordered Ordinal factor analysis. Not implemented yet; only `FALSE`
#'   is accepted.
#' @param sample_nobs,sample_mean,sample_cov,sample_cor Summary statistics
#'   used instead of `data`: the sample size, the item means (raw model only),
#'   and the covariance or correlation matrix. The covariance and raw models
#'   need `sample_cov`.
#' @param model Data model: `"cor"` (default) fits the correlation matrix of
#'   the standardized items, `"cov"` their covariance matrix, and `"raw"` the
#'   means and the covariance matrix.
#' @param lambda_prior Prior on the loadings: `"unit_vector"` (default), the
#'   unit-vector prior of Rey-Sáez et al. (2026), or `"normal"`, independent
#'   normal priors in the metric of the items (not available for `"cor"`).
#' @param missing Treatment of missing values: `"listwise"` (default) keeps
#'   the complete rows, and `"FIML"` uses every observed value. FIML needs
#'   `data`.
#' @param method Estimation method: `"mcmc"` (default) samples with NUTS,
#'   `"advi"` fits automatic differentiation variational inference, and
#'   `"pathfinder"` runs Pathfinder.
#' @param backend Stan interface: `"auto"` (default) takes the first one
#'   installed among stanli, cmdstanr, stanr and rstan that runs `method`.
#'   stanli needs no C++ compiler; the others compile each model once and keep
#'   it for later sessions (see [befa_compile_models()]).
#' @param rotate `"varimax"` (default) aligns the draws with Varimax and E-RSP,
#'   and orders the factors by the variance they explain (F1 explains the
#'   most); `"none"` returns the loadings as Stan sampled them.
#' @param rsp_args List with `max_iter` (default 1000) and `threshold`
#'   (default 1e-6) for [rsp_align()].
#' @param prior Named list of prior hyperparameters that replace the defaults
#'   (see Details): `xi`, `h2`, `lambda`, `psi`, `sigma` and `nu`.
#' @param factor_scores Logical. Draw the factor scores of every person (see
#'   [befa_factor_scores()]); needs `data`. The draws follow `seed`.
#' @param compute_fit_indices Logical. Compute the Bayesian fit indices and
#'   PSIS-LOO (see [befa_fit_measures()]); needs `data`. Outside the
#'   correlation model, this fits the null model with the same backend and
#'   method.
#' @param compute_reliability Logical. Compute the omega coefficients (see
#'   [befa_reliability()]).
#' @param loo_args List with `r_eff` and `cores` for the PSIS-LOO of the fit
#'   indices.
#' @param verbose Logical. Show the progress bar, the convergence checks and
#'   other messages.
#' @param ... Options of the estimation method. MCMC takes `chains` (default
#'   4), `parallel_chains` (chains that run at once; default up to 4, leaving
#'   two cores free), `iter_warmup` and `iter_sampling` (default 1000 each),
#'   `adapt_delta` (0.9), `max_treedepth` (15), `seed`, `refresh` (iterations
#'   between progress updates; `0` hides the progress bar) and `init`; the
#'   rstan names `iter` (warmup included), `warmup`, `cores` and `control` are
#'   also accepted. ADVI takes `algorithm` (`"meanfield"` or `"fullrank"`),
#'   `iter`, `draws`, `tol_rel_obj`, `seed`, `refresh` and `init`; Pathfinder
#'   takes `draws`, `num_paths`, `iter`, `psis_resample`, `seed`, `refresh`
#'   and `init`.
#'
#' @details
#' **_The factor model_**
#'
#' Each observed vector is \eqn{\mathbf{y}_i = \boldsymbol{\nu} +
#' \mathbf{\Lambda} \boldsymbol{\eta}_i + \boldsymbol{\epsilon}_i}, with
#' \eqn{\boldsymbol{\eta}_i \sim \mathcal{N}(\mathbf{0}, \mathbf{I}_M)} and
#' \eqn{\boldsymbol{\epsilon}_i \sim \mathcal{N}(\mathbf{0}, \mathbf{\Psi})},
#' so that \eqn{\mathbf{y}_i \sim \mathcal{N}(\boldsymbol{\nu},
#' \mathbf{\Lambda}\mathbf{\Lambda}^\top + \mathbf{\Psi})}. The correlation
#' model fits the standardized items (\eqn{\boldsymbol{\nu} = \mathbf{0}}, unit
#' variances), the covariance model the centered items, and the raw model
#' also estimates the means.
#'
#' **_Prior distributions_**
#'
#' * Unit-vector prior: each row of the standardized loadings is
#'   \eqn{\sqrt{h^2_j}\,\mathbf{z}_j / \lVert \mathbf{z}_j \rVert}, with a
#'   direction uniform on the sphere and a communality
#'   \eqn{h^2_j \sim \text{Beta}(a, b)} (`h2`, default `c(1, 1)`, which makes
#'   the joint prior on the loadings uniform). The norm of \eqn{\mathbf{z}_j}
#'   has density proportional to \eqn{\lVert \mathbf{z}_j \rVert^{\xi}
#'   e^{-\lVert \mathbf{z}_j \rVert^2/2}} (`xi`, default 100), which keeps it
#'   away from zero; it does not affect the loadings. With one factor, each
#'   loading follows a Beta prior scaled to (-1, 1) (`lambda`, default
#'   `c(1, 1)`). The uniquenesses are \eqn{1 - h^2_j}, so Heywood cases cannot
#'   occur. Outside the correlation model, the item standard deviations
#'   \eqn{\sigma_j} follow a half-Student-t (`sigma`, default `c(3, 0, 2.5)`).
#' * Normal prior (covariance and raw models): \eqn{\lambda_{jm} \sim
#'   \mathcal{N}(m s_j, s\, s_j)} and \eqn{\psi_j \sim \text{Inv-Gamma}(a,
#'   b\, s_j^2)}, where \eqn{s_j} is the standard deviation of item j (`lambda`,
#'   default `c(0, 10)`; `psi`, default `c(0.5, 0.5)`). Scaling by \eqn{s_j}
#'   makes the prior independent of the units of the items.
#' * Means (raw model): \eqn{\nu_j \sim \mathcal{N}(m, s)} (`nu`, default
#'   `c(0, 40)`).
#'
#' **_Estimation_**
#'
#' The Stan program is generated for the structure of the model and compiled
#' once per backend. Its likelihood uses the Woodbury identity, whose cost
#' grows with the number of items times the number of factors instead of the
#' cube of the number of items; under FIML, the rows are grouped by the items
#' they observe. Each chain runs in its own background process, starting from
#' the maximum-likelihood EFA that [stats::factanal()] computes, and the
#' convergence checks use the aligned parameters only.
#'
#' @return An object of class `befa` with:
#' * `$draws`: A `draws_array` of the final parameters in the metric of the
#'   model: `Lambda_std`, `Psi_std` and `h2` for the correlation model;
#'   `Lambda`, `Psi`, `sigma`, `h2` and (raw model) `nu` for the unit-vector
#'   prior; `Lambda`, `Psi` and (raw model) `nu` for the normal prior. Use
#'   [extract_posterior_draws()] with `metric = "std"` for the standardized
#'   solution.
#' * `$diagnostics`: The convergence checks (see the messages after the fit).
#' * `$fit`: The method, backend, settings, sampler diagnostics and run times,
#'   and the draws of the parameters of the Stan program with their log
#'   densities. Those parameters are not identified and are not meant for
#'   inference; [befa_bridge_sampler()] needs them.
#' * `$rsp`: The reference loadings, objective and sign and permutation
#'   vectors of the alignment, and `factor_order`, the order in which the
#'   aligned factors were placed (`NULL` with `rotate = "none"`).
#' * `$spec`, `$stan_code`, `$stan_data`, `$data`: The structure of the model,
#'   its Stan program, the data passed to Stan and the prepared data.
#' * `$eta`: The draws of the factor scores (if `factor_scores = TRUE`), also
#'   returned by [extract_posterior_draws()] in every metric.
#' * `$reliability`: The omega coefficients (if `compute_reliability = TRUE`).
#' * `$fit_indices`: The fit indices (if `compute_fit_indices = TRUE`).
#' * `$n_factors`, `$model_type`, `$lambda_prior`, `$missing`, `$method`,
#'   `$backend`, `$rotation`, `$priors_used`, `$options` and `$call`.
#'
#' @references
#' Rey-Sáez, R., Franco-Martínez, A., Revuelta, J., & Vadillo, M. A. (2026).
#' A unified framework for psychometrics in experimental psychology: The
#' standardized generalized hierarchical factor model. *Psychological
#' Methods*. Advance online publication. \doi{10.1037/met0000846}
#'
#' Rey-Sáez, R., & Revuelta, J. (2026). An efficient rotation-sign-permutation
#' algorithm to solve rotational indeterminacy in Bayesian exploratory factor
#' analysis. *Bayesian Analysis*. Advance publication.
#' \doi{10.1214/26-BA1611}
#'
#' @examples
#' \dontrun{
#' fit <- befa(HS_data, n_factors = 3, seed = 17)
#' extract_posterior_draws(fit, pars = "Lambda_std", format = "rvars")
#' posterior_summaries(fit, pars = "h2")
#' }
#'
#' @export
befa <- function(data = NULL, n_factors = NULL, ordered = FALSE,
                 sample_nobs = NULL, sample_mean = NULL, sample_cov = NULL,
                 sample_cor = NULL, model = c("cor", "cov", "raw"),
                 lambda_prior = c("unit_vector", "normal"),
                 missing = c("listwise", "FIML"),
                 method = c("mcmc", "advi", "pathfinder"), backend = "auto",
                 rotate = c("varimax", "none"), rsp_args = NULL, prior = list(),
                 factor_scores = FALSE, compute_fit_indices = TRUE,
                 compute_reliability = TRUE, loo_args = NULL, verbose = TRUE,
                 ...) {
  call <- match.call()
  model <- match.arg(model)
  lambda_prior <- match.arg(lambda_prior)
  missing <- match.arg(missing)
  method <- match.arg(method)
  rotate <- match.arg(rotate)
  .check_befa_args(data, n_factors, ordered, sample_nobs, sample_cov,
                   sample_cor, missing, factor_scores, compute_fit_indices,
                   compute_reliability, verbose)

  # Data, structure and Stan inputs
  prepared <- .prepare_befa_data(data, model, missing, sample_nobs, sample_cov,
                                 sample_cor, sample_mean)
  if (n_factors >= prepared$J) {
    .befa_stop("`n_factors` must be smaller than the number of items (",
               prepared$J, ").")
  }
  spec <- .befa_spec(model, lambda_prior, n_factors, .data_patterns(prepared),
                     .choose_engine())
  priors <- prepare_befa_priors(prior, model, lambda_prior, n_factors)
  stan_data <- .stan_data(spec, prepared, n_factors, priors)
  settings <- .sampler_settings(method, list(...))
  rsp_config <- validate_rsp_args(rsp_args)
  loo_config <- validate_loo_args(loo_args)
  backend <- .resolve_backend(backend, method, verbose = verbose)
  inits <- .befa_inits(settings, spec, prepared, n_factors, priors)

  # Estimation, alignment and convergence checks
  code <- .stan_code(spec)
  fit <- .run_stan(code, stan_data, inits, backend, settings, quiet = !verbose)
  post <- .befa_posterior(fit, spec, prepared$J, n_factors, rotate, rsp_config)
  diagnostics <- .befa_diagnostics(post$draws, fit, settings)
  if (verbose && settings$refresh > 0) .befa_report(diagnostics, settings)
  .befa_warn_issues(diagnostics)

  object <- structure(list(
    draws = post$draws, diagnostics = diagnostics, fit = fit, rsp = post$rsp,
    spec = spec, stan_code = code, stan_data = stan_data, data = prepared,
    n_factors = as.integer(n_factors), model_type = model,
    lambda_prior = lambda_prior, missing = missing, method = method,
    backend = backend, rotation = rotate, priors_used = priors,
    options = list(factor_scores = factor_scores, ordered = ordered,
                   rsp_args = rsp_config, loo_args = loo_config),
    call = call
  ), class = "befa")
  .befa_downstream(object, factor_scores, compute_reliability,
                   compute_fit_indices, settings$seed, verbose)
}

# Internal functions

# This function adds to a fit the quantities computed from its final draws:
# the factor scores (drawn from `seed`, so a fit is reproducible), the omega
# coefficients and the fit indices. The scores and the fit indices need the
# raw data; with summary statistics they are skipped with a warning. A fit
# index that fails leaves the fit intact and is reported as a warning.
.befa_downstream <- function(object, factor_scores, compute_reliability,
                             compute_fit_indices, seed, verbose) {
  summary_data <- object$data$summary
  if (factor_scores && summary_data) {
    .befa_warn("Factor scores need the raw `data`; they are skipped with ",
               "summary statistics.")
  } else if (factor_scores) {
    if (verbose) cli::cli_alert_info("Drawing the factor scores.")
    object$eta <- withr::with_seed(seed, .factor_score_draws(object))
  }
  if (compute_reliability) object$reliability <- befa_reliability(object)
  if (compute_fit_indices && summary_data) {
    .befa_warn("Fit indices need the raw `data`; they are skipped with ",
               "summary statistics.")
  } else if (compute_fit_indices) {
    if (verbose) {
      cli::cli_alert_info("Computing the fit indices and the null model.")
    }
    object$fit_indices <- tryCatch(befa_fit_measures(object), error = function(e) {
      .befa_warn("The fit indices could not be computed: ", conditionMessage(e))
      NULL
    })
  }
  object
}

# This function returns the starting values of a fit: those given in `init`
# (one list per chain or path, or one list for all), or the ones built from
# the ML EFA. MCMC needs one list per chain, Pathfinder one per path and ADVI
# a single list.
.befa_inits <- function(settings, spec, prepared, n_factors, priors) {
  n <- switch(settings$method, mcmc = settings$chains, advi = 1L,
              pathfinder = settings$num_paths)
  init <- settings$init
  if (is.null(init)) {
    init <- .start_values(spec, prepared, n_factors, priors, n, settings$seed)
  } else if (!is.list(init)) {
    .befa_stop("`init` must be a list of starting values, or one such list ",
               "per chain.")
  } else if (!is.null(names(init))) {
    init <- rep(list(init), n)
  } else if (length(init) != n) {
    .befa_stop("`init` has ", length(init), " lists of starting values, but ",
               "the fit needs ", n, ".")
  }
  if (settings$method == "advi") init[[1L]] else init
}
