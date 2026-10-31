/*
 * Pointwise log-likelihood of the null (independence) model
 *
 * Author: Ricardo Rey-Sáez
 * email: ricardoreysaez95@gmail.com
 * Modification date: 31/10/2026
 *
 * The null model has a diagonal covariance matrix, so the FIML log-density of
 * a person is the sum of univariate normal log-densities over the items that
 * person answered. This replaces the R loop over draws of v0.1.0.
 */

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <cmath>

using namespace arma;

static const double HALF_LOG_2PI = 0.5 * std::log(2.0 * arma::datum::pi);

// [[Rcpp::export]]
arma::mat null_loglik_cpp(const arma::mat& Y,           // N x J, NA = missing
                          const arma::mat& nu_draws,    // S x J item means
                          const arma::mat& sigma_draws  // S x J item SDs
) {
  const int N = Y.n_rows;
  const int J = Y.n_cols;
  const int S = sigma_draws.n_rows;

  if ((int)sigma_draws.n_cols != J || (int)nu_draws.n_cols != J ||
      (int)nu_draws.n_rows != S) {
    Rcpp::stop("nu_draws and sigma_draws must both be S x J, with J = %d.", J);
  }

  // Missing values arrive from R as NaN; they add nothing to the density
  const arma::umat observed = (Y == Y);

  arma::mat ll(S, N, arma::fill::zeros);
  for (int s = 0; s < S; ++s) {
    for (int j = 0; j < J; ++j) {
      const double mu = nu_draws(s, j);
      const double sd = sigma_draws(s, j);
      const double log_sd = std::log(sd);
      for (int i = 0; i < N; ++i) {
        if (!observed(i, j)) continue;
        const double z = (Y(i, j) - mu) / sd;
        ll(s, i) += -HALF_LOG_2PI - log_sd - 0.5 * z * z;
      }
    }
  }
  return ll;
}
