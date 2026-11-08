# befa-posterior.R
# Turn the Stan draws into aligned parameters in the metric of the model
# Author: Ricardo Rey-Sáez
# Last modified: 08-11-2026

# Internal functions

# This function returns the final draws of a fit and the result of the
# alignment. The standardized loadings are rebuilt from the Stan parameters,
# aligned with rsp_align() (unless `rotate = "none"`), and then returned in
# the metric of the model: standardized for the correlation model, in the
# units of the items for the covariance and raw models.
.befa_posterior <- function(fit, spec, J, M, rotate = "varimax",
                            rsp_args = validate_rsp_args(NULL)) {
  dims <- dim(fit$draws)
  n_iter <- dims[1L]
  n_chains <- dims[2L]
  std <- .standardized_parameters(fit$draws, spec, J, M)

  rsp <- NULL
  lambda <- std$lambda
  if (rotate != "none") {
    aligned <- rsp_align(lambda, n_items = J, n_factors = M, n_chains = n_chains,
                         format = "column_major", max_iter = rsp_args$max_iter,
                         threshold = rsp_args$threshold)
    lambda <- unname(aligned$Lambda_hat_mcmc)
    attr(lambda, "nchains") <- NULL
    rsp <- aligned[c("Lambda_star", "objective", "sign_vectors", "perm_vectors")]
    # The alignment leaves the factors in an arbitrary order. They are put in
    # the order of the variance they explain, as psych does, so F1 is the
    # strongest factor in every output.
    rsp$factor_order <- .factor_order(lambda, J, M)
    lambda <- lambda[, .factor_columns(rsp$factor_order, J), drop = FALSE]
    rsp$Lambda_star <- rsp$Lambda_star[, rsp$factor_order, drop = FALSE]
  }

  draws <- .model_metric(lambda, std, spec, J, M)
  list(draws = .as_chains(draws, n_iter, n_chains), rsp = rsp)
}

# This function returns the order of the factors by the variance of the
# standardized items they explain: the posterior mean of the sum of their
# squared loadings, largest first.
.factor_order <- function(lambda, J, M) {
  order(colSums(matrix(colMeans(lambda^2), J, M)), decreasing = TRUE)
}

# This function returns the columns of a column-major matrix of loadings
# draws (J per factor) that put its factors in the given order.
.factor_columns <- function(order, J) {
  unlist(lapply(order, function(m) (m - 1L) * J + seq_len(J)), use.names = FALSE)
}

# This function returns, from the draws of the Stan parameters, the
# standardized loadings (draws stacked chain after chain, one column per
# loading in column-major order), the standardized uniquenesses, the item
# standard deviations that undo the standardization (NULL in the correlation
# model) and the item means of the raw model.
.standardized_parameters <- function(draws, spec, J, M) {
  get <- function(par) .stacked_draws(draws, par)
  out <- list(scale = NULL, nu = NULL)
  if (spec$loadings == "uv") {
    Z <- get("Z")
    h2 <- get("h2")
    r2 <- .item_sums(Z^2, J, M)
    out$lambda <- Z * .repeat_items(sqrt(h2 / r2), M)
    out$psi <- 1 - h2
  } else if (spec$loadings == "uni") {
    out$lambda <- get("lambda_std")
    out$psi <- 1 - out$lambda^2
  } else {
    Lambda <- get("Lambda")
    psi <- get("psi")
    variance <- .item_sums(Lambda^2, J, M) + psi
    out$scale <- sqrt(variance)
    out$lambda <- Lambda / .repeat_items(out$scale, M)
    out$psi <- psi / variance
  }
  if (.spec_has_scale(spec)) out$scale <- get("sigma")
  if (spec$model_type == "raw") out$nu <- get("nu")
  out
}

# This function returns the final parameters in the metric of the model, with
# Stan-style names: Lambda_std, Psi_std and h2 in the correlation model;
# Lambda, Psi, sigma (unit-vector prior), h2 (unit-vector prior) and nu (raw
# model) otherwise.
.model_metric <- function(lambda, std, spec, J, M) {
  h2 <- 1 - std$psi
  if (spec$model_type == "cor") {
    parts <- list(Lambda_std = lambda, Psi_std = std$psi, h2 = h2)
  } else {
    parts <- list(Lambda = lambda * .repeat_items(std$scale, M),
                  Psi = std$psi * std$scale^2)
    if (spec$loadings != "normal") {
      parts$sigma <- std$scale
      parts$h2 <- h2
    }
    parts$nu <- std$nu
  }
  .named_draws(parts, J, M)
}

# This function converts the final draws of a fit to another metric:
# "model" returns them as they are, "std" standardizes the loadings and
# uniquenesses (and adds the communalities), and "raw" asks for the units of
# the items, which the correlation model does not have.
.metric_draws <- function(object, metric = c("model", "std", "raw")) {
  metric <- match.arg(metric)
  draws <- object$draws
  if (object$model_type == "cor" && metric == "raw") {
    .befa_stop("The correlation model has no metric of the items; use ",
               "`metric = \"std\"` or fit the covariance or raw model.")
  }
  if (metric == "std" && object$model_type != "cor") {
    J <- object$stan_data$J
    M <- object$n_factors
    Lambda <- .stacked_draws(draws, "Lambda")
    Psi <- .stacked_draws(draws, "Psi")
    variance <- .item_sums(Lambda^2, J, M) + Psi
    parts <- list(Lambda_std = Lambda / .repeat_items(sqrt(variance), M),
                  Psi_std = Psi / variance, h2 = 1 - Psi / variance)
    dims <- dim(draws)
    draws <- .as_chains(.named_draws(parts, J, M), dims[1L], dims[2L])
  }
  # The factor scores are standardized in every metric.
  if (!is.null(object$eta)) draws <- posterior::bind_draws(draws, object$eta)
  draws
}

# This function returns, draw by draw (chains stacked), the loadings,
# uniquenesses and means in the metric of the data the model was fitted to:
# standardized for the correlation model, in the units of the items
# otherwise, with zero means outside the raw model. The factor scores and the
# fit indices are computed from them.
.fitted_parameters <- function(object) {
  draws <- object$draws
  J <- object$stan_data$J
  std <- object$model_type == "cor"
  lambda <- .stacked_draws(draws, if (std) "Lambda_std" else "Lambda")
  psi <- .stacked_draws(draws, if (std) "Psi_std" else "Psi")
  nu <- if (object$model_type == "raw") {
    .stacked_draws(draws, "nu")
  } else {
    matrix(0, nrow(psi), J)
  }
  dims <- dim(draws)
  list(lambda = lambda, psi = psi, nu = nu, J = J, M = object$n_factors,
       n_iter = dims[1L], n_chains = dims[2L])
}

# This function returns the draws of one Stan parameter as a plain matrix,
# with the chains stacked one after another and the columns in Stan's
# (column-major) order.
.stacked_draws <- function(draws, par) {
  out <- posterior::as_draws_matrix(posterior::subset_draws(draws, variable = par))
  out <- unclass(out)
  attributes(out) <- list(dim = dim(out))
  out
}

# This function sums, for every item, the M columns (one per factor) that a
# column-major matrix of loadings holds for it.
.item_sums <- function(x, J, M) {
  Reduce(`+`, lapply(seq_len(M), function(m) x[, (m - 1L) * J + seq_len(J), drop = FALSE]))
}

# This function repeats a matrix with one column per item M times, so it
# lines up with a column-major matrix of loadings.
.repeat_items <- function(x, M) {
  do.call(cbind, rep(list(x), M))
}

# This function joins draws of named parameters into one matrix with Stan
# names: matrices of loadings become name[j,m] and vectors name[j].
.named_draws <- function(parts, J, M) {
  parts <- parts[!vapply(parts, is.null, TRUE)]
  blocks <- lapply(names(parts), function(par) {
    x <- parts[[par]]
    colnames(x) <- if (ncol(x) == J * M && grepl("^Lambda", par)) {
      paste0(par, "[", rep(seq_len(J), M), ",", rep(seq_len(M), each = J), "]")
    } else {
      paste0(par, "[", seq_len(J), "]")
    }
    x
  })
  do.call(cbind, blocks)
}

# This function turns draws stacked chain after chain back into a draws_array
# with their chains.
.as_chains <- function(x, n_iter, n_chains) {
  posterior::as_draws_array(array(x, dim = c(n_iter, n_chains, ncol(x)),
                                  dimnames = list(NULL, NULL, colnames(x))))
}
