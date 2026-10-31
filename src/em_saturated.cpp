/*
 * EM estimation of the saturated-model moments with missing data
 *
 * Author: Ricardo Rey-Sáez
 * email: ricardoreysaez95@gmail.com
 * Modification date: 31/10/2026
 *
 * Same algorithm as the R version of v0.1.0, itself based on lavaan's
 * lav_mvnorm_missing_h1_estimate_moments(): start from the observed column
 * means and a diagonal covariance, then alternate the E-step (expected
 * sufficient statistics given the observed values) and the M-step until the
 * largest change falls below tol. The C++ version groups the rows by their
 * pattern of missing values, so each pattern needs one matrix inverse per
 * iteration instead of one per row.
 */

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <map>
#include <vector>

using namespace arma;

// Rows that share one pattern of observed items
struct MissingPattern {
  arma::uvec obs;
  arma::uvec mis;
  arma::uvec rows;
};

static std::vector<MissingPattern> missing_patterns(const arma::umat& observed) {
  std::map<std::vector<int>, std::vector<arma::uword>> groups;
  for (arma::uword i = 0; i < observed.n_rows; ++i) {
    std::vector<int> key(observed.n_cols);
    for (arma::uword j = 0; j < observed.n_cols; ++j) key[j] = observed(i, j);
    groups[key].push_back(i);
  }

  std::vector<MissingPattern> out;
  for (const auto& g : groups) {
    MissingPattern p;
    arma::uvec flag = arma::conv_to<arma::uvec>::from(g.first);
    p.obs = arma::find(flag == 1);
    p.mis = arma::find(flag == 0);
    p.rows = arma::conv_to<arma::uvec>::from(g.second);
    out.push_back(p);
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List em_h1_moments_cpp(const arma::mat& Y, int max_iter = 500,
                             double tol = 1e-5) {
  const arma::uword N = Y.n_rows;
  const arma::uword J = Y.n_cols;

  // Missing values arrive from R as NaN
  const arma::umat observed = (Y == Y);
  const std::vector<MissingPattern> patterns = missing_patterns(observed);

  // Start: observed column means and a diagonal covariance
  arma::vec Mu(J);
  arma::vec var0(J);
  for (arma::uword j = 0; j < J; ++j) {
    arma::vec yj = Y.col(j);
    arma::vec yo = yj.elem(arma::find(observed.col(j)));
    // An item nobody answered starts at mean 0 and variance 1 (v0.1.0 gave NaN)
    Mu(j) = yo.n_elem > 0 ? arma::mean(yo) : 0.0;
    var0(j) = yo.n_elem > 0 ? arma::mean(arma::square(yo - Mu(j))) : 1.0;
    if (!std::isfinite(var0(j)) || var0(j) == 0.0) var0(j) = 1.0;
  }
  arma::mat Sigma = arma::diagmat(var0);

  int iterations = 0;
  for (int iter = 0; iter < max_iter; ++iter) {
    iterations = iter + 1;
    const arma::vec Mu0 = Mu;
    const arma::mat Sigma0 = Sigma;

    // E-step: expected sums (T1) and cross-products (T2)
    arma::vec T1(J, arma::fill::zeros);
    arma::mat T2(J, J, arma::fill::zeros);

    for (const MissingPattern& p : patterns) {
      // Rows with no observed item add nothing, as in v0.1.0
      if (p.obs.n_elem == 0) continue;

      arma::mat Yp = Y.rows(p.rows);

      if (p.mis.n_elem == 0) {
        T1 += arma::sum(Yp, 0).t();
        T2 += Yp.t() * Yp;
        continue;
      }

      // Regression of the missing items on the observed ones
      const arma::mat S_mo = Sigma.submat(p.mis, p.obs);
      const arma::mat B = S_mo * arma::inv(Sigma.submat(p.obs, p.obs));
      const arma::mat cond_cov = Sigma.submat(p.mis, p.mis) - B * S_mo.t();

      for (arma::uword r = 0; r < p.rows.n_elem; ++r) {
        arma::vec y_full = Yp.row(r).t();
        arma::vec y_obs = y_full.elem(p.obs);
        arma::vec y_complete(J);
        y_complete.elem(p.obs) = y_obs;
        y_complete.elem(p.mis) = Mu.elem(p.mis) + B * (y_obs - Mu.elem(p.obs));

        T1 += y_complete;
        T2 += y_complete * y_complete.t();
        T2.submat(p.mis, p.mis) += cond_cov;
      }
    }

    // M-step
    Mu = T1 / static_cast<double>(N);
    Sigma = T2 / static_cast<double>(N) - Mu * Mu.t();

    // Keep Sigma positive definite with a tiny ridge, as in v0.1.0
    arma::vec ev = arma::eig_sym(Sigma);
    if (arma::any(ev < 1e-6)) {
      Sigma.diag() += Sigma.diag().max() * 1e-8;
    }

    const double delta = std::max(arma::abs(Mu - Mu0).max(),
                                  arma::abs(Sigma - Sigma0).max());
    if (delta < tol) break;
  }

  return Rcpp::List::create(
    Rcpp::Named("Mu") = Rcpp::NumericVector(Mu.begin(), Mu.end()),
    Rcpp::Named("Sigma") = Sigma,
    Rcpp::Named("iterations") = iterations
  );
}
