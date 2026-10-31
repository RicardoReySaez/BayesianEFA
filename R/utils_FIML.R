#' Compute FIML Log-Likelihood (R equivalent of Stan's incomplete_obs_lpdf)
#'
#' For each observation, computes the MVN log-density using only observed
#' variables (marginalizing over missing ones). This is the proper FIML approach.
#' It is highly optimized to compute complete cases vectorized.
#'
#' @param data_fiml A list containing `Y_comp` (matrix of complete cases),
#'   `Y_miss` (matrix of incomplete cases), `comp_idx` (original indices for complete cases),
#'   and `miss_idx` (original indices for incomplete cases).
#' @param mu Vector. Mean vector (length J).
#' @param Sigma Matrix. Covariance matrix (J x J).
#' @return Numeric vector of log-likelihoods matching the original N.
#' @keywords internal
#' @noRd
compute_fiml_loglik <- function(data_fiml, mu, Sigma) {
  # Total number of observations is the sum of complete and incomplete cases
  N_comp <- length(data_fiml$comp_idx)
  N_miss <- length(data_fiml$miss_idx)
  N_total <- N_comp + N_miss

  ll_vec <- numeric(N_total)

  # 1. Fast-Path for complete cases (Vectorized C++ call)
  if (N_comp > 0) {
    ll_vec[data_fiml$comp_idx] <- mvnfast::dmvn(data_fiml$Y_comp, mu = mu, sigma = Sigma, log = TRUE)
  }

  # 2. Row-by-row FIML for incomplete cases
  if (N_miss > 0) {
    Y_miss <- data_fiml$Y_miss

    # Pre-allocate subset vector for missing cases
    ll_miss <- numeric(N_miss)

    for (i in seq_len(N_miss)) {
      obs_idx <- which(!is.na(Y_miss[i, ]))
      n_obs <- length(obs_idx)

      if (n_obs == 0) {
        ll_miss[i] <- 0 # No observed data -> no contribution
      } else {
        # Marginalize over observed dimensions
        y_obs <- Y_miss[i, obs_idx]
        mu_obs <- mu[obs_idx]
        Sigma_obs <- Sigma[obs_idx, obs_idx, drop = FALSE]

        ll_miss[i] <- mvnfast::dmvn(matrix(y_obs, nrow = 1), mu = mu_obs, sigma = Sigma_obs, log = TRUE)
      }
    }
    # Assign the evaluated missing case likelihoods back to their original indices
    ll_vec[data_fiml$miss_idx] <- ll_miss
  }

  return(ll_vec)
}


# EM algorithm for the saturated-model moments with missing data

# This function estimates the means and covariance matrix of the saturated
# model from data with missing values, which keeps the chi-square comparisons
# valid under FIML. It runs the EM algorithm of lavaan
# (lav_mvnorm_missing_h1_estimate_moments) in C++, grouping the rows by their
# pattern of missing values.
estimate_h1_moments_em <- function(Y, max_iter = 500L, tol = 1e-05) {
  res <- em_h1_moments_cpp(unname(as.matrix(Y)), as.integer(max_iter), tol)
  list(Mu = res$Mu, Sigma = res$Sigma)
}
