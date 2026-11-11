# bench-precision.R
# Accuracy of the log density of each engine as a uniqueness goes to zero
# Author: Ricardo Rey-Sáez
# Last modified: 31-10-2026

# The Woodbury identity subtracts two terms of order 1/psi, so it could lose
# digits in a near-Heywood case. This benchmark evaluates the log density and
# its gradient of both engines with stanli at the same point, with the
# smallest uniqueness going from 1e-1 to 1e-12, and compares them with a
# reference whose log-likelihood is computed with 120-bit MPFR arithmetic.
# It uses the normal prior, whose uniquenesses are Stan parameters with a
# lower bound: their log transform round-trips to a relative 1e-16, so the
# point Stan evaluates is the point of the reference.

source("bench/bench-helpers.R")
source("tests/testthat/helper-stan.R")
out <- "bench/results/precision.csv"
unlink(out)
bits <- 120

# This function returns the Cholesky factor of an MPFR matrix.
.chol_mpfr <- function(A) {
  n <- nrow(A)
  L <- Rmpfr::mpfrArray(0, bits, dim = c(n, n))
  for (j in seq_len(n)) {
    s <- A[j, j]
    if (j > 1) for (k in seq_len(j - 1)) s <- s - L[j, k]^2
    L[j, j] <- sqrt(s)
    if (j < n) for (i in (j + 1):n) {
      s <- A[i, j]
      if (j > 1) for (k in seq_len(j - 1)) s <- s - L[i, k] * L[j, k]
      L[i, j] <- s / L[j, j]
    }
  }
  L
}

# This function returns the log-likelihood of the prepared (centered) data in
# MPFR arithmetic, pattern by pattern: each pattern's Cholesky factor gives
# its determinant, and forward substitution its quadratic forms.
.ll_mpfr <- function(prepared, Lambda, psi) {
  Lm <- Rmpfr::mpfr(Lambda, bits)
  total <- Rmpfr::mpfr(0, bits)
  for (g in prepared$groups) {
    o <- g$obs
    Sigma <- Lm[o, , drop = FALSE] %*% t(Lm[o, , drop = FALSE])
    for (j in seq_along(o)) Sigma[j, j] <- Sigma[j, j] + Rmpfr::mpfr(psi[o[j]], bits)
    L <- .chol_mpfr(Sigma)
    X <- Rmpfr::mpfr(t(prepared$Y[g$rows, o, drop = FALSE]), bits)
    Z <- X
    for (i in seq_along(o)) {
      s <- X[i, ]
      if (i > 1) for (k in seq_len(i - 1)) s <- s - L[i, k] * Z[k, ]
      Z[i, ] <- s / L[i, i]
    }
    n <- length(g$rows)
    logdet <- 2 * Reduce(`+`, lapply(seq_along(o), function(i) log(L[i, i])))
    total <- total - 0.5 * (n * length(o) * log(2 * Rmpfr::Const("pi", bits)) +
                              n * logdet + sum(Z^2))
  }
  total
}

# This function returns the log-likelihood of the prepared data in double
# precision, row by row, as .reference_lp() computes it.
.ll_double <- function(prepared, Sigma) {
  sum(vapply(seq_len(nrow(prepared$Y)), function(i) {
    o <- which(!is.na(prepared$Y[i, ]))
    mvnfast::dmvn(prepared$Y[i, o], rep(0, length(o)), Sigma[o, o, drop = FALSE],
                  log = TRUE)
  }, 0))
}

J <- 10
M <- 2
for (missing in c(0, 0.1)) {
  Y <- .bench_data(J, M, N = 500, missing = missing, seed = 3) * 1.5
  prepared <- .prepare_befa_data(Y, "cov", if (missing) "FIML" else "listwise")
  prior <- prepare_befa_priors(list(), "cov", "normal", M, J)
  item_sd <- .item_sd(prepared)
  set.seed(4)
  Lambda <- matrix(stats::runif(J * M, 0.2, 0.7), J, M) * item_sd
  for (psi_min in 10^-(1:12)) {
    psi <- (stats::runif(J, 0.3, 0.6) * item_sd^2)
    psi[1] <- psi_min
    pars <- list(Lambda = Lambda, psi = psi)
    Sigma <- tcrossprod(Lambda) + diag(psi)
    spec_c <- .befa_spec("cov", "normal", M, .data_patterns(prepared), "cholesky")
    reference <- .reference_lp(spec_c, prepared, pars, prior, item_sd) -
      .ll_double(prepared, Sigma) +
      Rmpfr::asNumeric(.ll_mpfr(prepared, Lambda, psi))
    res <- list()
    for (engine in c("cholesky", "woodbury")) {
      spec <- .befa_spec("cov", "normal", M, .data_patterns(prepared), engine)
      m <- stanli::stanli_model(code = .stan_code(spec),
                                data = .stan_data(spec, prepared, M, prior))
      res[[engine]] <- stanli::log_prob_grad(m, stanli::unconstrain(m, pars))
    }
    g_c <- res$cholesky$grad
    g_w <- res$woodbury$grad
    row <- data.frame(
      patterns = length(prepared$groups), psi_min = psi_min,
      lp_reference = reference,
      error_cholesky = res$cholesky$lp - reference,
      error_woodbury = res$woodbury$lp - reference,
      grad_rel_diff = max(abs(g_w - g_c)) / max(abs(g_c))
    )
    .bench_write(row, out)
    cat(sprintf("P=%3d psi_min=%.0e  chol %+.2e  wood %+.2e  grad %.1e\n",
                row$patterns, psi_min, row$error_cholesky, row$error_woodbury,
                row$grad_rel_diff))
  }
}
