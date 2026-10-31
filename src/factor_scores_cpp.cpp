/*
 * Bayesian Factor Score computation via RcppArmadillo
 *
 * Author: Ricardo Rey-Sáez
 * email: ricardoreysaez95@gmail.com
 * Modification date: 31/10/2026
 *
 * Draws the factor scores of every person from their conditional posterior
 * given each draw of the loadings, uniquenesses and means, in the metric of
 * the data the model was fitted to. Incomplete rows condition on their
 * observed items only (FIML). The arithmetic and the order of the random
 * numbers are those of v0.1.0, which read Sigma instead of the uniquenesses.
 */

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::export]]
arma::mat compute_factor_scores_cpp(
    const arma::mat& lambda,            // (S x J*M) loadings, column-major per draw
    const arma::mat& psi,               // (S x J) uniquenesses
    const arma::mat& nu,                // (S x J) item means (zeros if none)
    const arma::mat& data,              // (N x J) data in the metric of the model
    int J, int M,
    bool has_missing,
    const Rcpp::List& obs_patterns      // list of integer vectors (0-based observed indices per person)
) {

  const int S = lambda.n_rows;
  const int N = data.n_rows;
  const int P = N * M;

  // Output: one row per draw, persons within factors (column-major eta[i,m])
  arma::mat eta_out(S, P, arma::fill::zeros);

  const arma::mat I_M = arma::eye<arma::mat>(M, M);

  for (int s = 0; s < S; ++s) {

    // ── Parameters of this draw ──
    arma::mat L_s = arma::reshape(lambda.row(s).t(), J, M);
    arma::vec Psi_inv = 1.0 / psi.row(s).t();
    arma::vec nu_s = nu.row(s).t();

    // ── Complete-data path: vectorized over all N persons ──
    if (!has_missing) {
      // L_s' * Psi^{-1} * L_s  (M x M)
      arma::mat L_weighted = L_s.each_col() % Psi_inv;  // (J x M) with rows scaled
      arma::mat LtPsiInvL  = L_s.t() * L_weighted;      // (M x M)

      // Sigma_eta = (I_M + L' Psi^{-1} L)^{-1}
      arma::mat Sigma_eta = arma::inv_sympd(I_M + LtPsiInvL);

      // Cholesky of Sigma_eta for sampling
      arma::mat chol_Sig;
      bool chol_ok = arma::chol(chol_Sig, Sigma_eta);
      if (!chol_ok) {
        chol_Sig = arma::chol(Sigma_eta);
      }

      // Conditional means for all persons: (N x M)
      arma::mat Y_centered = data.each_row() - nu_s.t();
      arma::mat mu_all = Y_centered * L_weighted * Sigma_eta;

      // Sample: eta = mu + Z * chol(Sigma_eta)
      arma::mat Z = arma::randn<arma::mat>(N, M);
      eta_out.row(s) = arma::vectorise(mu_all + Z * chol_Sig).t();

    // ── Missing-data path: person by person ──
    } else {
      for (int i = 0; i < N; ++i) {
        Rcpp::IntegerVector obs_i_R = obs_patterns[i];
        int n_obs_i = obs_i_R.size();

        if (n_obs_i == 0) {
          // No observed data: sample from prior N(0, I)
          arma::vec eta_prior = arma::randn<arma::vec>(M);
          for (int m = 0; m < M; ++m) eta_out(s, m * N + i) = eta_prior(m);
          continue;
        }

        arma::uvec obs_i(n_obs_i);
        for (int k = 0; k < n_obs_i; ++k) {
          obs_i(k) = static_cast<arma::uword>(obs_i_R[k]);
        }

        // Loadings, means and precisions of the observed items
        arma::mat L_i       = L_s.rows(obs_i);
        arma::vec nu_i      = nu_s.elem(obs_i);
        arma::vec psi_inv_i = Psi_inv.elem(obs_i);
        arma::vec y_full    = data.row(i).t();
        arma::vec y_i       = y_full.elem(obs_i);

        // Conditional posterior
        arma::mat L_w_i       = L_i.each_col() % psi_inv_i;
        arma::mat Sigma_eta_i = arma::inv_sympd(I_M + L_i.t() * L_w_i);
        arma::vec mu_i        = Sigma_eta_i * (L_w_i.t() * (y_i - nu_i));

        arma::mat chol_i;
        bool chol_ok = arma::chol(chol_i, Sigma_eta_i);
        if (!chol_ok) {
          chol_i = arma::chol(Sigma_eta_i + 1e-8 * I_M);
        }
        arma::vec eta_i = mu_i + chol_i.t() * arma::randn<arma::vec>(M);

        for (int m = 0; m < M; ++m) eta_out(s, m * N + i) = eta_i(m);
      }
    }
  }

  return eta_out;
}
