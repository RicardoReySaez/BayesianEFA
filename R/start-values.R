# start-values.R
# Starting values of the chains from a maximum-likelihood EFA
# Author: Ricardo Rey-Sáez
# Last modified: 31-10-2026

# The chains start where stats::factanal(), lavaan and psych::fa(fm = "ml")
# start their own estimation: from the maximum-likelihood EFA of the
# correlation matrix, whose uniquenesses begin at (1 - M / (2J)) / diag(R^-1).

# Internal functions

# This function returns the maximum-likelihood EFA of a correlation matrix with
# the algorithm of stats::factanal(): the uniquenesses are optimized with
# L-BFGS-B between 0.005 and 1, and the loadings follow from them in closed
# form, through the eigen-decomposition of Psi^-1/2 R Psi^-1/2. It is written
# here because factanal.fit.mle() is not exported from stats.
.ml_efa <- function(R, M) {
  J <- ncol(R)
  # This function returns the loadings implied by the uniquenesses psi.
  loadings <- function(psi) {
    sc <- 1 / sqrt(psi)
    E <- eigen(R * tcrossprod(sc), symmetric = TRUE)
    L <- E$vectors[, seq_len(M), drop = FALSE] %*%
      diag(sqrt(pmax(E$values[seq_len(M)] - 1, 0)), M)
    sqrt(psi) * L
  }
  # This function returns the ML discrepancy as factanal() writes it.
  objective <- function(psi) {
    sc <- 1 / sqrt(psi)
    e <- eigen(R * tcrossprod(sc), symmetric = TRUE,
               only.values = TRUE)$values[-seq_len(M)]
    -(sum(log(e) - e) - M + J)
  }
  # This function returns the gradient of the discrepancy.
  gradient <- function(psi) {
    L <- loadings(psi)
    diag(tcrossprod(L) + diag(psi, J) - R) / psi^2
  }

  start <- (1 - 0.5 * M / J) / diag(solve(R))
  fit <- stats::optim(start, objective, gradient, method = "L-BFGS-B",
                      lower = 0.005, upper = 1,
                      control = list(fnscale = 1, parscale = rep(0.01, J)))
  list(loadings = loadings(fit$par), uniquenesses = fit$par,
       converged = fit$convergence == 0L)
}

# This function returns a starting solution in the standardized metric: the
# ML EFA, or, when it fails, the one-step principal-axis solution that psych
# starts from, with the squared multiple correlations as communalities.
.start_solution <- function(R, M) {
  efa <- tryCatch(.ml_efa(R, M), error = function(e) NULL)
  if (!is.null(efa) && efa$converged) return(efa)
  smc <- 1 - 1 / diag(solve(R))
  reduced <- R
  diag(reduced) <- smc
  E <- eigen(reduced, symmetric = TRUE)
  L <- E$vectors[, seq_len(M), drop = FALSE] %*%
    diag(sqrt(pmax(E$values[seq_len(M)], 0)), M)
  list(loadings = L, uniquenesses = 1 - pmin(rowSums(L^2), 0.995),
       converged = FALSE)
}

# This function returns the correlation matrix, standard deviations and
# means of the prepared data: from the summary statistics, from the complete
# rows, or from the EM saturated model when values are missing.
.start_moments <- function(prepared) {
  if (prepared$summary) {
    S <- prepared$S
    return(list(R = stats::cov2cor(S), sd = sqrt(diag(S)),
                mean = if (is.null(prepared$y_bar)) rep(0, prepared$J) else prepared$y_bar))
  }
  Y <- prepared$Y
  if (anyNA(Y)) {
    h1 <- estimate_h1_moments_em(Y)
    return(list(R = stats::cov2cor(h1$Sigma), sd = sqrt(diag(h1$Sigma)),
                mean = h1$Mu))
  }
  n <- nrow(Y)
  S <- stats::cov(Y) * (n - 1) / n
  list(R = stats::cov2cor(S), sd = sqrt(diag(S)), mean = colMeans(Y))
}

# This function returns one list of starting values per chain, on the
# constrained scale of the Stan parameters of `spec`. Every chain starts from
# the same EFA solution with its own small perturbation, drawn from
# `seed + chain` so that the user's random stream is left untouched. The
# unit-vector rows of Z start at the mode of their norm, sqrt(xi + M - 1),
# where the repulsion places them.
.start_values <- function(spec, prepared, n_factors, prior, chains, seed) {
  moments <- .start_moments(prepared)
  J <- prepared$J
  M <- as.integer(n_factors)
  sol <- if (spec$loadings != "none") .start_solution(moments$R, M)
  h2 <- if (!is.null(sol)) pmin(pmax(rowSums(sol$loadings^2), 0.05), 0.95)
  s <- moments$sd

  lapply(seq_len(chains), function(chain) {
    withr::with_seed(seed + chain, {
      init <- list()
      if (spec$loadings == "uv") {
        dir <- sol$loadings + matrix(stats::rnorm(J * M, 0, 0.05), J, M)
        radius <- sqrt(prior$xi + M - 1)
        init$h2 <- pmin(pmax(h2 + stats::rnorm(J, 0, 0.02), 0.05), 0.95)
        init$Z <- radius * dir / sqrt(rowSums(dir^2))
      } else if (spec$loadings == "uni") {
        init$lambda_std <- pmin(pmax(sol$loadings[, 1L] + stats::rnorm(J, 0, 0.02),
                                     -0.95), 0.95)
      } else if (spec$loadings == "normal") {
        init$Lambda <- s * (sol$loadings + matrix(stats::rnorm(J * M, 0, 0.05), J, M))
        init$psi <- s^2 * (1 - h2) * exp(stats::rnorm(J, 0, 0.05))
      }
      if (.spec_has_scale(spec)) init$sigma <- s * exp(stats::rnorm(J, 0, 0.05))
      if (spec$model_type == "raw") {
        init$nu <- moments$mean + stats::rnorm(J, 0, 0.05) * s
      }
      init
    })
  })
}
