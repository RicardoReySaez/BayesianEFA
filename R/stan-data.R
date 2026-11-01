# stan-data.R
# Prepare the data of a BEFA model and translate them for the Stan code
# Author: Ricardo Rey-Sáez
# Last modified: 31-10-2026

# Internal functions

# This function prepares the observations the model is fitted to. Raw data
# are standardized (correlation model) or centered (covariance model) with
# the moments of the observed values, as in v0.1.0, and their rows are grouped
# by the items they observe. Listwise deletion keeps the complete rows only;
# under FIML, rows without any observed value add nothing and are dropped.
# Summary statistics give a single group of `sample_nobs` rows, described by
# their maximum-likelihood covariance (or correlation) matrix and means.
.prepare_befa_data <- function(data = NULL, model_type = c("cor", "cov", "raw"),
                               missing = c("listwise", "FIML"),
                               sample_nobs = NULL, sample_cov = NULL,
                               sample_cor = NULL, sample_mean = NULL) {
  model_type <- match.arg(model_type)
  missing <- match.arg(missing)
  if (is.null(data)) {
    return(.prepare_summary_data(model_type, sample_nobs, sample_cov,
                                 sample_cor, sample_mean))
  }

  Y <- unname(as.matrix(data))
  .check_matrix(Y, "data")
  N_in <- nrow(Y)
  keep <- if (missing == "listwise") {
    stats::complete.cases(Y)
  } else {
    rowSums(!is.na(Y)) > 0L
  }
  if (!all(keep)) {
    why <- if (missing == "listwise") "with missing values" else "without any observed value"
    .befa_warn(sum(!keep), " of ", N_in, " rows ", why, " were removed.")
  }
  Y <- Y[keep, , drop = FALSE]
  if (nrow(Y) == 0L) .befa_stop("No rows are left to fit the model.")

  # Standardize or center each item with its observed values
  centers <- colMeans(Y, na.rm = TRUE)
  if (model_type != "raw") Y <- sweep(Y, 2L, centers)
  if (model_type == "cor") Y <- sweep(Y, 2L, apply(Y, 2L, .ml_sd), "/")

  list(Y = Y, N = nrow(Y), J = ncol(Y), groups = .missing_patterns(Y),
       summary = FALSE, model_type = model_type)
}

# This function prepares summary statistics: the sample size, the
# maximum-likelihood covariance matrix (or the correlation matrix in the
# correlation model) and, in the raw model, the means.
.prepare_summary_data <- function(model_type, sample_nobs, sample_cov,
                                  sample_cor, sample_mean) {
  if (is.null(sample_nobs)) {
    .befa_stop("Summary statistics need `sample_nobs`.")
  }
  .check_count(sample_nobs, "sample_nobs")
  if (model_type != "cor" && is.null(sample_cov)) {
    .befa_stop("The ", model_type, " model needs `sample_cov`; a correlation ",
               "matrix carries no variances.")
  }
  if (model_type == "raw" && is.null(sample_mean)) {
    .befa_stop("The raw model needs `sample_mean` with summary statistics.")
  }

  use_cor <- model_type == "cor" && !is.null(sample_cor)
  S <- if (use_cor) sample_cor else sample_cov
  .check_matrix(S, if (use_cor) "sample_cor" else "sample_cov")
  N <- as.integer(sample_nobs)
  S <- unname(if (model_type == "cor") stats::cov2cor(S) else S * (N - 1) / N)
  J <- ncol(S)
  list(Y = NULL, N = N, J = J,
       groups = list(list(rows = seq_len(N), obs = seq_len(J))),
       summary = TRUE, S = S,
       y_bar = if (model_type == "raw") as.numeric(sample_mean),
       model_type = model_type)
}

# This function tells whether every row observes every item, which the Stan
# code treats as a single pattern.
.data_patterns <- function(prepared) {
  one <- length(prepared$groups) == 1L &&
    length(prepared$groups[[1L]]$obs) == prepared$J
  if (one) "single" else "multiple"
}

# This function groups the rows of Y by the items they observe. The group of
# complete rows comes first and the others follow in a fixed order, so the
# same data always give the same Stan data.
.missing_patterns <- function(Y) {
  observed <- !is.na(Y)
  key <- apply(observed, 1L, function(r) paste(as.integer(r), collapse = ""))
  levels <- unique(key)
  complete <- paste(rep("1", ncol(Y)), collapse = "")
  levels <- c(intersect(complete, levels), sort(setdiff(levels, complete)))
  lapply(levels, function(k) {
    rows <- which(key == k)
    list(rows = rows, obs = which(observed[rows[1L], ]))
  })
}

# This function returns the maximum-likelihood standard deviation of the
# observed values of x, which makes the cross-products of standardized data
# equal N times the sample correlation matrix.
.ml_sd <- function(x) {
  x <- x[!is.na(x)]
  sqrt(mean((x - mean(x))^2))
}

# This function returns a root R of a symmetric positive semidefinite matrix
# S (R %*% t(R) = S): its Cholesky factor or, when S is singular, the scaled
# eigenvectors of its positive eigenvalues. A null matrix gives one column of
# zeros, since the Stan code needs at least one column. A Cholesky factor
# with a pivot near zero also counts as singular: whether chol() fails on a
# singular matrix or ends with a pivot of rounding size depends on the
# platform's LAPACK.
.matrix_root <- function(S) {
  R <- tryCatch(t(chol(S)), error = function(e) NULL)
  if (!is.null(R) && min(diag(R))^2 > max(diag(S), 0) * 1e-12) return(R)
  e <- eigen(S, symmetric = TRUE)
  keep <- e$values > max(e$values, 0) * 1e-12
  if (!any(keep)) return(matrix(0, nrow(S), 1L))
  e$vectors[, keep, drop = FALSE] %*% diag(sqrt(e$values[keep]), sum(keep))
}

# This function returns a root of the scatter matrix t(Y) %*% Y of the rows
# of Y with as few columns as possible: the rows themselves when there are
# fewer rows than items, and the root of the scatter matrix otherwise.
.scatter_root <- function(Y) {
  if (nrow(Y) < ncol(Y)) return(t(Y))
  .matrix_root(crossprod(Y))
}

# This function returns the standard deviation of each item, which scales the
# normal prior on the loadings: the maximum-likelihood one with complete data
# or summary statistics, and that of the EM saturated model with missing
# values.
.item_sd <- function(prepared) {
  if (prepared$summary) return(sqrt(diag(prepared$S)))
  if (anyNA(prepared$Y)) {
    return(sqrt(diag(estimate_h1_moments_em(prepared$Y)$Sigma)))
  }
  apply(prepared$Y, 2L, .ml_sd)
}

# This function returns the data list that the Stan code of `spec` reads:
# dimensions, the root of the scatter matrix of each pattern (around its
# means in the raw model), the means of each pattern in the raw model, and
# the prior hyperparameters.
.stan_data <- function(spec, prepared, n_factors, prior) {
  if (spec$patterns != .data_patterns(prepared)) {
    .befa_stop("The model structure expects ", spec$patterns, " patterns of ",
               "observed items, but the data have the other kind.")
  }
  raw <- spec$model_type == "raw"
  out <- list(N = as.integer(prepared$N), J = as.integer(prepared$J))
  if (spec$loadings %in% c("uv", "normal")) out$M <- as.integer(n_factors)

  # One scatter root (and mean vector) per pattern
  pieces <- lapply(prepared$groups, function(g) {
    if (prepared$summary) {
      return(list(n = prepared$N, obs = g$obs, y_bar = prepared$y_bar,
                  root = .matrix_root(prepared$N * prepared$S)))
    }
    Yg <- prepared$Y[g$rows, g$obs, drop = FALSE]
    y_bar <- if (raw) colMeans(Yg)
    if (raw) Yg <- sweep(Yg, 2L, y_bar)
    list(n = length(g$rows), obs = g$obs, y_bar = y_bar,
         root = .scatter_root(Yg))
  })

  if (spec$patterns == "single") {
    out$K <- ncol(pieces[[1L]]$root)
    out$S_root <- pieces[[1L]]$root
    if (raw) out$y_bar <- pieces[[1L]]$y_bar
  } else {
    roots <- lapply(pieces, `[[`, "root")
    out$P <- length(pieces)
    out$n_p <- vapply(pieces, `[[`, 1L, "n")
    out$J_p <- vapply(pieces, function(x) length(x$obs), 1L)
    out$K_p <- vapply(roots, ncol, 1L)
    out$obs_idx <- as.integer(unlist(lapply(pieces, `[[`, "obs")))
    out$n_root <- as.integer(sum(lengths(roots)))
    out$root_flat <- unlist(lapply(roots, as.vector))
    if (raw) out$ybar_flat <- unlist(lapply(pieces, `[[`, "y_bar"))
  }

  item_sd <- if (spec$loadings == "normal") .item_sd(prepared)
  c(out, .stan_prior_data(spec, prior, item_sd))
}

# This function translates the resolved prior into the hyperparameters that
# the Stan code reads. Under the normal prior, the loadings of item j are
# scaled by its standard deviation s_j and the scale of its inverse-gamma
# uniqueness by s_j^2 (if psi ~ IG(a, b), then c psi ~ IG(a, c b)), so the
# prior describes the standardized solution and does not depend on the units
# of the items.
.stan_prior_data <- function(spec, prior, item_sd = NULL) {
  out <- list()
  if (spec$loadings == "uv") {
    out$pr_h2 <- prior$h2
    out$pr_xi <- prior$xi
  }
  if (spec$loadings == "uni") out$pr_lambda <- prior$lambda
  if (spec$loadings == "normal") {
    out$pr_lambda_loc <- prior$lambda[1L] * item_sd
    out$pr_lambda_scale <- prior$lambda[2L] * item_sd
    out$pr_psi_shape <- prior$psi[1L]
    out$pr_psi_scale <- prior$psi[2L] * item_sd^2
  }
  if (.spec_has_scale(spec)) out$pr_sigma <- prior$sigma
  if (spec$model_type == "raw") out$pr_nu <- prior$nu
  out
}
