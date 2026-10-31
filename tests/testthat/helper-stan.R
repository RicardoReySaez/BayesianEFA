# helper-stan.R
# Shared helpers for the tests of the generated Stan code
# Author: Ricardo Rey-Sáez
# Last modified: 31-10-2026

# This function returns every model structure the generator supports, named
# by its label: the correlation, covariance and raw models with each prior on
# the loadings and each kind of pattern, and the null models.
.all_specs <- function() .spec_grid()

# This function returns the number of factors a structure is tested with.
.spec_factors <- function(spec) if (spec$loadings == "uni") 1L else 2L

# This function simulates raw data from a two-factor model, with item means
# and SDs away from 0 and 1. With `missing = TRUE` some values are removed,
# including a row with a single observed item, so several patterns appear and
# some of them have fewer rows than items.
.simulated_data <- function(N = 40L, J = 5L, missing = FALSE, seed = 41) {
  set.seed(seed)
  Lambda <- matrix(c(0.8, 0.7, 0.6, 0.1, 0.2, 0.1, 0.2, 0.1, 0.7, 0.8), J, 2)
  Sigma <- tcrossprod(Lambda) + diag(0.4, J)
  sds <- seq(0.8, 2, length.out = J)
  Y <- mvnfast::rmvn(N, mu = seq(1, 3, length.out = J),
                     sigma = Sigma * tcrossprod(sds))
  if (missing) {
    Y[sample(N * J, round(0.12 * N * J))] <- NA
    Y[1L, ] <- c(NA, NA, Y[1L, 3L], NA, NA)
    Y[2L, ] <- c(Y[2L, 1L], NA, NA, NA, NA)
    Y[rowSums(!is.na(Y)) == 0L, 1L] <- 1
  }
  Y
}

# This function prepares the data a structure is tested with: complete data
# for a single pattern, data with missing values (FIML) otherwise.
.prepared_for <- function(spec, seed = 41) {
  multiple <- spec$patterns == "multiple"
  .prepare_befa_data(.simulated_data(missing = multiple, seed = seed),
                     spec$model_type, if (multiple) "FIML" else "listwise")
}

# This function returns the default prior of a structure.
.prior_for <- function(spec) {
  lambda_prior <- if (spec$loadings == "normal") "normal" else "unit_vector"
  n_factors <- if (spec$loadings == "uni") 1L else 2L
  prepare_befa_priors(list(), spec$model_type, lambda_prior, n_factors)
}

# This function draws values of the parameters of a structure, on their
# constrained scale.
.random_pars <- function(spec, J, M, seed = 42) {
  set.seed(seed)
  pars <- switch(spec$loadings,
    uv = list(h2 = stats::runif(J, 0.2, 0.8),
              Z = matrix(stats::rnorm(J * M, 0, 3), J, M)),
    uni = list(lambda_std = stats::runif(J, -0.8, 0.8)),
    normal = list(Lambda = matrix(stats::rnorm(J * M, 0, 0.8), J, M),
                  psi = stats::runif(J, 0.3, 1.2)),
    none = list()
  )
  if (.spec_has_scale(spec)) pars$sigma <- stats::runif(J, 0.6, 1.8)
  if (spec$model_type == "raw") pars$nu <- stats::rnorm(J, 2, 0.5)
  pars
}

# This function returns the normalizing constant of the density
# ||z||^xi exp(-||z||^2 / 2) on R^M.
.repulsion_log_const <- function(xi, M) {
  0.5 * M * log(pi) + 0.5 * (xi + M) * log(2) + lgamma(0.5 * (xi + M)) -
    lgamma(0.5 * M)
}

# This function returns the log posterior density, with the Jacobian of the
# constrained parameters, that the Stan code of a structure must give. It is
# computed independently of the Stan code: the likelihood comes row by row
# from the prepared data (or from the summary statistics), not from scatter
# roots.
.reference_lp <- function(spec, prepared, pars, prior, item_sd = NULL) {
  J <- prepared$J
  lp <- 0

  if (spec$loadings == "uv") {
    r2 <- rowSums(pars$Z^2)
    L <- pars$Z * sqrt(pars$h2 / r2)
    Sigma <- tcrossprod(L) + diag(1 - pars$h2)
    lp <- lp + sum(stats::dbeta(pars$h2, prior$h2[1], prior$h2[2], log = TRUE)) +
      sum(0.5 * prior$xi * log(r2) - 0.5 * r2) -
      J * .repulsion_log_const(prior$xi, ncol(pars$Z)) +
      sum(log(pars$h2) + log(1 - pars$h2))
  } else if (spec$loadings == "uni") {
    l <- pars$lambda_std
    Sigma <- tcrossprod(l) + diag(1 - l^2)
    lp <- lp + sum(stats::dbeta((l + 1) / 2, prior$lambda[1], prior$lambda[2],
                                log = TRUE) - log(2)) +
      sum(log((1 + l) * (1 - l) / 2))
  } else if (spec$loadings == "normal") {
    Sigma <- tcrossprod(pars$Lambda) + diag(pars$psi)
    a <- prior$psi[1]
    b <- prior$psi[2] * item_sd^2
    lp <- lp +
      sum(stats::dnorm(pars$Lambda, prior$lambda[1] * item_sd,
                       prior$lambda[2] * item_sd, log = TRUE)) +
      sum(a * log(b) - lgamma(a) - (a + 1) * log(pars$psi) - b / pars$psi) +
      sum(log(pars$psi))
  } else {
    Sigma <- diag(J)
  }

  if (.spec_has_scale(spec)) {
    s <- pars$sigma
    df <- prior$sigma[1]
    loc <- prior$sigma[2]
    sc <- prior$sigma[3]
    Sigma <- Sigma * tcrossprod(s)
    lp <- lp + sum(stats::dt((s - loc) / sc, df, log = TRUE) - log(sc)) -
      J * stats::pt((0 - loc) / sc, df, lower.tail = FALSE, log.p = TRUE) +
      sum(log(s))
  }

  mu <- rep(0, J)
  if (spec$model_type == "raw") {
    mu <- pars$nu
    lp <- lp + sum(stats::dnorm(mu, prior$nu[1], prior$nu[2], log = TRUE))
  }

  if (prepared$summary) {
    N <- prepared$N
    d <- (if (is.null(prepared$y_bar)) rep(0, J) else prepared$y_bar) - mu
    ll <- -0.5 * N * (J * log(2 * pi) + determinant(Sigma)$modulus[1] +
                        sum(diag(solve(Sigma, prepared$S))) +
                        sum(d * solve(Sigma, d)))
  } else {
    ll <- sum(vapply(seq_len(nrow(prepared$Y)), function(i) {
      o <- which(!is.na(prepared$Y[i, ]))
      mvnfast::dmvn(prepared$Y[i, o], mu[o], Sigma[o, o, drop = FALSE],
                    log = TRUE)
    }, 0))
  }
  lp + ll
}
