/*
 * Varimax rotation of posterior draws of the loading matrix
 *
 * Author: Ricardo Rey-Sáez
 * email: ricardoreysaez95@gmail.com
 * Modification date: 31/10/2026
 *
 * Rotates every draw with the same algorithm, stopping rule and defaults as
 * stats::varimax(x, normalize = FALSE, eps = 1e-5), so rsp_align() keeps the
 * results described in the paper while it no longer loops over draws in R.
 */

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace arma;

// Varimax rotation of one J x M loading matrix (no Kaiser normalization).
// It follows stats::varimax line by line: the criterion is the sum of the
// singular values of B, and the loop stops once it grows by less than a
// relative eps, or after max_iter iterations.
static arma::mat varimax_one(const arma::mat& L, const double eps,
                             const int max_iter) {
  const int p  = L.n_rows;
  const int nc = L.n_cols;

  arma::mat TT = arma::eye<arma::mat>(nc, nc);
  arma::mat z, B, U, V;
  arma::vec sv;
  double d = 0.0;

  for (int i = 0; i < max_iter; ++i) {
    z = L * TT;
    // z %*% diag(colSums(z^2)) / p, as in stats::varimax
    arma::rowvec colsq = arma::sum(arma::square(z), 0);
    arma::mat zd = z.each_row() % colsq;
    B = L.t() * (arma::pow(z, 3) - zd / p);
    // La.svd() and the default svd() of Armadillo both call LAPACK's dgesdd
    arma::svd(U, sv, V, B);
    TT = U * V.t();
    const double dpast = d;
    d = arma::sum(sv);
    if (d < dpast * (1.0 + eps)) break;
  }
  return L * TT;
}

// [[Rcpp::export]]
arma::mat varimax_draws_cpp(const arma::mat& draws, int J, int M,
                            double eps = 1e-5, int max_iter = 1000) {
  if ((int)draws.n_cols != J * M) {
    Rcpp::stop("Each draw must have J * M = %d values, not %d.", J * M,
               (int)draws.n_cols);
  }

  // With a single factor there is nothing to rotate
  if (M < 2) return draws;

  const int n_draws = draws.n_rows;
  arma::mat out(n_draws, J * M);

  for (int s = 0; s < n_draws; ++s) {
    // Each row holds a column-major (factor by factor) loading matrix
    arma::mat L = arma::reshape(draws.row(s).t(), J, M);
    out.row(s) = arma::vectorise(varimax_one(L, eps, max_iter)).t();
  }
  return out;
}
