# helper-befa.R
# A befa object built without Stan, for the functions that read its draws
# Author: Ricardo Rey-Sáez
# Last modified: 11-11-2026

# This function returns a befa object whose standardized loadings scatter
# around the true loadings of .simulated_data(), in the metric of `model`:
# two chains of `n_iter` draws, with the prepared data, structure, priors,
# settings and diagnostics that befa() would store. No Stan program is run.
# The standardized draws do not depend on `model`, so fits of the same seed
# share them.
.fake_befa <- function(model = "cor", lambda_prior = "unit_vector",
                       missing = "listwise", N = 60L, n_iter = 100L, seed = 7,
                       M = 2L) {
  data <- .simulated_data(N = N, missing = missing == "FIML", seed = seed)
  prepared <- suppressWarnings(.prepare_befa_data(data, model, missing))
  J <- prepared$J
  S <- 2L * n_iter
  spec <- .befa_spec(model, lambda_prior, M, .data_patterns(prepared))
  priors <- prepare_befa_priors(list(), model, lambda_prior, M, J)

  set.seed(seed)
  true <- c(0.8, 0.7, 0.6, 0.1, 0.2, 0.1, 0.2, 0.1, 0.7, 0.8,
            0.05, -0.05, 0.05, -0.05, 0.05)[seq_len(J * M)]
  lambda <- matrix(true, S, J * M, byrow = TRUE) +
    matrix(stats::rnorm(S * J * M, 0, 0.03), S)
  std <- list(lambda = lambda, psi = 1 - .item_sums(lambda^2, J, M))
  if (model != "cor") {
    sds <- apply(prepared$Y, 2L, stats::sd, na.rm = TRUE)
    std$scale <- matrix(sds, S, J, byrow = TRUE) *
      exp(matrix(stats::rnorm(S * J, 0, 0.02), S))
  }
  if (model == "raw") {
    std$nu <- matrix(colMeans(prepared$Y, na.rm = TRUE), S, J, byrow = TRUE) +
      matrix(stats::rnorm(S * J, 0, 0.05), S)
  }
  draws <- .as_chains(.model_metric(lambda, std, spec, J, M), n_iter, 2L)

  settings <- .sampler_settings("mcmc", list(chains = 2, iter_warmup = 10,
                                             iter_sampling = n_iter, seed = seed))
  fit <- list(method = "mcmc", backend = "stanli", settings = settings,
              sampler = NULL)
  # The draws are independent, so posterior warns that their ESS is capped.
  diagnostics <- suppressWarnings(.befa_diagnostics(draws, fit, settings))
  structure(list(
    draws = draws, diagnostics = diagnostics,
    fit = fit, rsp = NULL, spec = spec, stan_code = NULL,
    stan_data = .stan_data(spec, prepared, M, priors), data = prepared,
    n_factors = M, model_type = model, lambda_prior = lambda_prior,
    missing = missing, method = "mcmc", backend = "stanli",
    rotation = "varimax", priors_used = priors,
    options = list(factor_scores = FALSE, ordered = FALSE, prior_only = FALSE,
                   rsp_args = validate_rsp_args(NULL),
                   loo_args = validate_loo_args(NULL)),
    call = quote(befa())
  ), class = "befa")
}

# This function returns fixed null-model draws in place of a Stan fit: two
# chains of `n_iter` draws of sigma, and of nu for the raw model.
.null_draws <- function(J, raw, n_iter = 5L, seed = 21) {
  set.seed(seed)
  vars <- c(paste0("sigma[", seq_len(J), "]"),
            if (raw) paste0("nu[", seq_len(J), "]"))
  values <- c(stats::runif(n_iter * 2 * J, 0.8, 1.6),
              if (raw) stats::rnorm(n_iter * 2 * J, 0, 0.3))
  posterior::as_draws_array(array(values, dim = c(n_iter, 2, length(vars)),
                                  dimnames = list(NULL, NULL, vars)))
}
