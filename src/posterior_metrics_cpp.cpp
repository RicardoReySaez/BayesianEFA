/*
 * C++ backend for compute_posterior_metrics (RcppArmadillo)
 *
 * Author: Ricardo Rey-Sáez
 * email: ricardoreysaez95@gmail.com
 * Modification date: 23/02/2026
 *
 * Computes pointwise log-likelihood, chi-square, chi-square replicated (PPP),
 * and SRMR for each MCMC posterior draw. This replaces the R-level loop in
 * utils_fitmeasures.R with a high-performance C++ implementation.
 */

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace arma;

// Pre-computed constant: log(2 * pi)
static const double LOG_2PI = std::log(2.0 * arma::datum::pi);

// ──────────────────────────────────────────────────────────────────────
// Helper functions
// ──────────────────────────────────────────────────────────────────────

// MVN log-density for a single observation (J-vector, cholesky approach)
static inline double mvn_logpdf(const arma::vec &y, const arma::vec &mu,
                                const arma::mat &L_chol, double log_det,
                                double half_J_log2pi) {
  arma::vec z = arma::solve(arma::trimatl(L_chol), y - mu, arma::solve_opts::fast);
  return -half_J_log2pi - log_det - 0.5 * arma::dot(z, z);
}

// Helper: MVN log-density for a matrix of observations (N x J)
static arma::vec mvn_logpdf_mat(const arma::mat &Y, const arma::vec &mu,
                                  const arma::mat &Sigma,
                                  double half_J_log2pi) {
  int N = Y.n_rows;
  arma::mat L_chol = arma::chol(Sigma, "lower");
  double log_det = arma::sum(arma::log(L_chol.diag()));

  arma::vec ll(N);
  for (int i = 0; i < N; ++i) {
    arma::vec yi = Y.row(i).t();
    ll(i) = mvn_logpdf(yi, mu, L_chol, log_det, half_J_log2pi);
  }
  return ll;
}


// FIML log-density for a single observation with missing data
static double fiml_logpdf_single(const arma::rowvec &yi_row,
                                 const arma::vec &mu, const arma::mat &Sigma,
                                 const arma::uvec &obs_i) {
  int n_obs = obs_i.n_elem;
  if (n_obs == 0)
    return 0.0;

  arma::vec y_obs(n_obs);
  arma::vec mu_obs(n_obs);
  arma::mat Sigma_obs(n_obs, n_obs);

  for (int a = 0; a < n_obs; ++a) {
    y_obs(a) = yi_row(obs_i(a));
    mu_obs(a) = mu(obs_i(a));
    for (int b = 0; b < n_obs; ++b) {
      Sigma_obs(a, b) = Sigma(obs_i(a), obs_i(b));
    }
  }

  double half_n_log2pi = 0.5 * n_obs * LOG_2PI;
  arma::mat L_chol = arma::chol(Sigma_obs, "lower");
  double log_det = arma::sum(arma::log(L_chol.diag()));
  arma::vec z = arma::solve(arma::trimatl(L_chol), y_obs - mu_obs,
                            arma::solve_opts::fast);
  return -half_n_log2pi - log_det - 0.5 * arma::dot(z, z);
}


// Compute FIML pointwise log-likelihood for all N observations
static arma::vec
compute_fiml_ll(const arma::mat &Y_comp, const arma::mat &Y_miss,
                const arma::uvec &comp_idx, // 0-based
                const arma::uvec &miss_idx, // 0-based
                const arma::vec &mu, const arma::mat &Sigma, int N_total,
                const Rcpp::List &obs_patterns_miss, double half_J_log2pi) {

  arma::vec ll(N_total, arma::fill::zeros);

  // 1. Fast-path: vectorized for complete cases
  int N_comp = comp_idx.n_elem;
  if (N_comp > 0) {
    arma::vec ll_comp = mvn_logpdf_mat(Y_comp, mu, Sigma, half_J_log2pi);
    for (int i = 0; i < N_comp; ++i) {
      ll(comp_idx(i)) = ll_comp(i);
    }
  }

  // 2. Row-by-row FIML for incomplete cases
  int N_miss = miss_idx.n_elem;
  if (N_miss > 0) {
    for (int i = 0; i < N_miss; ++i) {
      Rcpp::IntegerVector obs_i_R = obs_patterns_miss[i];
      arma::uvec obs_i(obs_i_R.size());
      for (int k = 0; k < obs_i_R.size(); ++k) {
        obs_i(k) = static_cast<arma::uword>(obs_i_R[k]);
      }
      ll(miss_idx(i)) = fiml_logpdf_single(Y_miss.row(i), mu, Sigma, obs_i);
    }
  }

  return ll;
}

// Compute sample correlation matrix from a matrix Y (N x J)
static arma::mat compute_cor(const arma::mat &Y) {
  int N = Y.n_rows;
  arma::vec means = arma::mean(Y, 0).t();
  arma::mat centered = Y.each_row() - means.t();
  arma::mat cov_mat = (centered.t() * centered) / (N - 1);
  arma::vec sd_inv = 1.0 / arma::sqrt(cov_mat.diag());
  arma::mat cor_mat = cov_mat % (sd_inv * sd_inv.t());
  return cor_mat;
}

// Compute sample covariance matrix with ML denominator (N, not N-1)
static arma::mat compute_cov_ml(const arma::mat &Y) {
  int N = Y.n_rows;
  arma::vec means = arma::mean(Y, 0).t();
  arma::mat centered = Y.each_row() - means.t();
  return (centered.t() * centered) / N;
}

// cov2cor (same as R's stats::cov2cor)
static arma::mat cov2cor(const arma::mat &S) {
  arma::vec d_inv = 1.0 / arma::sqrt(S.diag());
  return S % (d_inv * d_inv.t());
}

// Compute SRMR from model-implied and sample correlation matrices
static double compute_srmr(const arma::mat &Sample_Cor,
                           const arma::mat &Imp_Cor, int J) {
  double ss = 0.0;
  int count = 0;
  for (int j = 0; j < J; ++j) {
    for (int i = j; i < J; ++i) {
      double diff = Sample_Cor(i, j) - Imp_Cor(i, j);
      ss += diff * diff;
      ++count;
    }
  }
  return std::sqrt(ss / count);
}

// ──────────────────────────────────────────────────────────────────────
// Compute posterior metrics main function
// ──────────────────────────────────────────────────────────────────────
// [[Rcpp::export]]
Rcpp::List compute_posterior_metrics_cpp(
    const arma::cube &post_draws, // (n_iter x n_chains x n_vars)
    const arma::mat &data,        // (N x J) — already scaled for cor model
    const arma::uvec &sig_idx,    // 0-based indices for Sigma/Rho
    const arma::uvec &nu_idx,     // 0-based indices for nu (empty if cor/cov)
    const arma::vec &M_sat,       // saturated mean vector (J)
    const arma::mat &S_sat,       // saturated cov/cor matrix (J x J)
    double ll_saturated,          // pre-computed saturated LL
    int model_type,               // 1 = raw, 2 = cov, 3 = cor
    bool has_missing,
    const arma::mat &Y_comp,      // complete cases subset (N_comp x J)
    const arma::mat &Y_miss,      // incomplete cases subset (N_miss x J)
    const arma::uvec &comp_idx_0, // 0-based indices of complete cases in original data
    const arma::uvec &miss_idx_0, // 0-based indices of incomplete cases
    const Rcpp::List &obs_patterns_miss // list of 0-based observed indices for each miss row
) {

  const int n_iter = post_draws.n_rows;
  const int n_chains = post_draws.n_cols;
  const int n_draws = n_iter * n_chains;
  const int N = data.n_rows;
  const int J = data.n_cols;
  const bool has_nu = (nu_idx.n_elem > 0);

  // Pre-compute MVN normalization constant (constant across all draws)
  const double half_J_log2pi = 0.5 * J * LOG_2PI;

  // Pre-compute Sample_Cor for SRMR (constant across draws)
  arma::mat Sample_Cor;
  if (model_type == 1 || model_type == 2) { // raw or cov
    Sample_Cor = cov2cor(S_sat);
  } else { // cor
    Sample_Cor = S_sat;
  }

  // Output containers
  arma::cube log_lik_cube(n_iter, n_chains, N, arma::fill::zeros);
  arma::vec chisq_vec(n_draws, arma::fill::zeros);
  arma::vec chisq_rep_vec(n_draws, arma::fill::zeros);
  arma::vec srmr_vec(n_draws, arma::fill::zeros);

  // Main loop: chains x iterations
  int draw_idx = 0;
  for (int c = 0; c < n_chains; ++c) {
    for (int s = 0; s < n_iter; ++s) {
      draw_idx = c * n_iter + s;

      // Extract Sigma_s (J x J) from draws
      arma::vec sig_vec(sig_idx.n_elem);
      for (arma::uword k = 0; k < sig_idx.n_elem; ++k) {
        sig_vec(k) = post_draws(s, c, sig_idx(k));
      }
      arma::mat Sigma_s = arma::reshape(sig_vec, J, J);

      // Extract Nu_s (J) from draws
      arma::vec Nu_s(J, arma::fill::zeros);
      if (has_nu) {
        for (arma::uword k = 0; k < nu_idx.n_elem; ++k) {
          Nu_s(k) = post_draws(s, c, nu_idx(k));
        }
      }

      // Pointwise log-likelihood
      arma::vec ll_s = compute_fiml_ll(Y_comp, Y_miss, comp_idx_0, miss_idx_0, Nu_s, Sigma_s,
                                       N, obs_patterns_miss, half_J_log2pi);
      log_lik_cube.tube(s, c) = ll_s;

      // Chi-square: 2 * (LL_Saturated - LL_Model)
      chisq_vec(draw_idx) = 2.0 * (ll_saturated - arma::accu(ll_s));

      // Posterior Predictive P-value (replicate data)
      // Generate Y_rep ~ MVN(Nu_s, Sigma_s)
      arma::mat L_model = arma::chol(Sigma_s, "lower");
      arma::mat Z = arma::randn<arma::mat>(N, J);
      arma::mat Y_rep = Z * L_model.t();
      Y_rep.each_row() += Nu_s.t();

      // Log-lik of Y_rep under the model
      double ll_rep_model = arma::accu(mvn_logpdf_mat(Y_rep, Nu_s, Sigma_s, half_J_log2pi));

      // Saturated moments of Y_rep
      arma::vec M_sat_rep;
      arma::mat S_sat_rep;

      if (model_type == 3) { // cor
        S_sat_rep = compute_cor(Y_rep);
        M_sat_rep = arma::zeros<arma::vec>(J);
      } else if (model_type == 2) { // cov
        S_sat_rep = compute_cov_ml(Y_rep);
        M_sat_rep = arma::zeros<arma::vec>(J);
      } else { // raw
        S_sat_rep = compute_cov_ml(Y_rep);
        M_sat_rep = arma::mean(Y_rep, 0).t();
      }

      // Log-lik of Y_rep under saturated moments
      double ll_sat_rep = arma::accu(mvn_logpdf_mat(Y_rep, M_sat_rep, S_sat_rep, half_J_log2pi));
      chisq_rep_vec(draw_idx) = 2.0 * (ll_sat_rep - ll_rep_model);

      // SRMR
      arma::mat Imp_Cor;
      if (model_type == 1 || model_type == 2) { // raw or cov
        Imp_Cor = cov2cor(Sigma_s);
      } else { // cor
        Imp_Cor = Sigma_s;
      }
      srmr_vec(draw_idx) = compute_srmr(Sample_Cor, Imp_Cor, J);
    }
  }

  return Rcpp::List::create(Rcpp::Named("log_lik") = Rcpp::wrap(log_lik_cube),
                            Rcpp::Named("chisq") = Rcpp::wrap(chisq_vec),
                            Rcpp::Named("chisq_rep") = Rcpp::wrap(chisq_rep_vec),
                            Rcpp::Named("srmr") = Rcpp::wrap(srmr_vec));
}
