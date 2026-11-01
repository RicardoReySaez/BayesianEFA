// log_density.cpp
// Log density of the BEFA models in unconstrained coordinates, and the map
// from the coordinates of bridge sampling to those coordinates
// Author: Ricardo Rey-Sáez
// Last modified: 01-11-2026

// The density is the one the generated Stan programs compute, with every
// additive constant and the Jacobian of the constrained parameters, so that
// it integrates to the marginal likelihood. The likelihood uses the Woodbury
// identity in the vectorized form of the Stan code: every pattern of
// observed items shares one matrix product, and each pattern needs only the
// Cholesky factor of its M x M capacitance matrix. Complete data and summary
// statistics are the case of a single pattern.

#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <vector>

namespace {

const double LOG_2PI = std::log(2.0 * M_PI);
const double NEG_INF = -std::numeric_limits<double>::infinity();

// The structure, priors and data of a model.
struct Density {
  int loadings;            // 0 unit-vector, 1 one factor (scaled Beta), 2 normal
  bool scale, raw, prior_only;
  int J, M;
  // Priors
  double h2_a, h2_b, xi, log_C;          // unit-vector
  double lam_a, lam_b;                   // one factor
  arma::vec lam_loc, lam_scale;          // normal: loadings
  double psi_shape;
  arma::vec psi_scale;                   // normal: uniquenesses
  double sig_df, sig_loc, sig_scale;     // item standard deviations
  double nu_loc, nu_scale;               // item means
  // Data of the patterns
  arma::mat obs_mat;    // P x J, 1 for the items each pattern observes
  arma::mat root_pad;   // K_tot x J, columns of the scatter roots as rows
  arma::vec n_obs;      // rows that observe each item
  arma::vec r2_item;    // squared norm of the root entries of each item
  arma::mat ybar_pad;   // P x J, pattern means (raw model)
  arma::vec n_p;        // rows per pattern
  arma::uvec K_p;       // columns of each root
  double ll_const;
};

// The map from the coordinates of bridge sampling to the unconstrained
// coordinates of the model: the canonical loadings T (Z[ord, ] = T Q) in
// polar coordinates (log norm of each row and a chart of its direction) or
// in Cartesian ones (log of the diagonal and the free entries).
struct BridgeMap {
  int type;                  // 0 polar, 1 Cartesian
  int p, M, n_u;
  arma::uvec inv_ord;        // row j of Z is row inv_ord[j] of T
  arma::uvec load_idx;       // positions of the loadings in u (0-based)
  arma::uvec other_idx;      // positions of the other parameters
  std::vector<arma::vec> centers;  // chart centres of the rows past M
  std::vector<arma::mat> bases;    // orthonormal bases of their tangents
  double log_vol;
};

double log_inv_logit(double x) {
  return x < 0 ? x - std::log1p(std::exp(x)) : -std::log1p(std::exp(-x));
}

double log_sech(double x) {
  x = std::fabs(x);
  return std::log(2.0) - x - std::log1p(std::exp(-2.0 * x));
}

// log |d omega / d z| of a direction chart with k tangent coordinates.
double chart_log_jac(double zeta, double theta, int k, double theta_max) {
  double ratio = zeta < 1e-10 ? theta_max / 2.0 : std::sin(theta) / zeta;
  return (k - 1) * std::log(ratio) + std::log(theta_max / 2.0) +
    2.0 * log_sech(zeta / 2.0);
}

// The log-likelihood with the vectorized Woodbury identity, for
// Sigma = A A' + diag(u) and item means nu.
double woodbury_loglik(const Density& d, const arma::mat& A,
                       const arma::vec& u, const arma::vec& nu) {
  const int R = A.n_cols;
  const arma::uword P = d.n_p.n_elem;
  arma::vec inv_u = 1.0 / u;
  arma::mat B = A.each_col() % arma::sqrt(inv_u);
  arma::mat W = A.each_col() % inv_u;
  arma::mat BB(d.J, R * R);
  for (int k = 0; k < R; ++k) {
    for (int l = 0; l < R; ++l) BB.col(k * R + l) = B.col(k) % B.col(l);
  }
  arma::mat C_flat = d.obs_mat * BB;
  arma::mat Q = d.root_pad * W;
  arma::mat D, Q_d;
  arma::vec d2;
  if (d.raw) {
    D = d.ybar_pad - (d.obs_mat.each_row() % nu.t());
    Q_d = D * W;
    d2 = arma::square(D) * inv_u;
  }

  double logdet_C = 0.0, quad_C = 0.0, mean_C = 0.0;
  arma::uword pos = 0;
  arma::mat C(R, R), L_C;
  for (arma::uword p = 0; p < P; ++p) {
    C = arma::reshape(C_flat.row(p).t(), R, R);
    C.diag() += 1.0;
    if (!arma::chol(L_C, C, "lower")) return NEG_INF;
    logdet_C += d.n_p(p) * arma::accu(arma::log(L_C.diag()));
    arma::mat X = arma::solve(arma::trimatl(L_C),
                              Q.rows(pos, pos + d.K_p(p) - 1).t(),
                              arma::solve_opts::fast);
    quad_C += arma::accu(arma::square(X));
    if (d.raw) {
      arma::vec y = arma::solve(arma::trimatl(L_C), Q_d.row(p).t(),
                                arma::solve_opts::fast);
      mean_C += d.n_p(p) * arma::accu(arma::square(y));
    }
    pos += d.K_p(p);
  }
  double ll = d.ll_const + 0.5 * arma::dot(d.n_obs, arma::log(inv_u)) -
    logdet_C - 0.5 * (arma::dot(d.r2_item, inv_u) - quad_C);
  if (d.raw) ll -= 0.5 * (arma::dot(d.n_p, d2) - mean_C);
  return ll;
}

// The log density at unconstrained coordinates u, in the order of the
// parameters of the Stan program.
double log_density(const Density& d, const arma::vec& u) {
  const int J = d.J, M = d.M;
  arma::uword pos = 0;
  double lp = 0.0;
  arma::mat A;
  arma::vec un;

  if (d.loadings == 0) {
    arma::vec x = u.subvec(pos, pos + J - 1);
    pos += J;
    arma::mat Z = arma::reshape(u.subvec(pos, pos + J * M - 1), J, M);
    pos += J * M;
    arma::vec log_h2(J), log_1m(J);
    for (int j = 0; j < J; ++j) {
      log_h2(j) = log_inv_logit(x(j));
      log_1m(j) = log_inv_logit(-x(j));
    }
    arma::vec h2 = arma::exp(log_h2);
    arma::vec r2 = arma::sum(arma::square(Z), 1);
    A = Z.each_col() % arma::sqrt(h2 / r2);
    un = arma::exp(log_1m);
    lp += arma::accu((d.h2_a - 1.0) * log_h2 + (d.h2_b - 1.0) * log_1m) -
      J * R::lbeta(d.h2_a, d.h2_b);
    lp += arma::accu(0.5 * d.xi * arma::log(r2) - 0.5 * r2) - J * d.log_C;
    lp += arma::accu(log_h2 + log_1m);
  } else if (d.loadings == 1) {
    arma::vec x = u.subvec(pos, pos + J - 1);
    pos += J;
    arma::vec lam(J), jac(J);
    for (int j = 0; j < J; ++j) {
      double lp1 = log_inv_logit(x(j)), lp0 = log_inv_logit(-x(j));
      lam(j) = std::exp(lp1) - std::exp(lp0);
      jac(j) = std::log(2.0) + lp1 + lp0;
    }
    A = lam;
    un = 1.0 - arma::square(lam);
    lp += arma::accu((d.lam_a - 1.0) * arma::log1p(lam) +
                     (d.lam_b - 1.0) * arma::log1p(-lam)) -
      J * (d.lam_a + d.lam_b - 1.0) * std::log(2.0) -
      J * R::lbeta(d.lam_a, d.lam_b);
    lp += arma::accu(jac);
  } else {
    arma::mat Lam = arma::reshape(u.subvec(pos, pos + J * M - 1), J, M);
    pos += J * M;
    arma::vec x = u.subvec(pos, pos + J - 1);
    pos += J;
    arma::vec psi = arma::exp(x);
    A = Lam;
    un = psi;
    for (int m = 0; m < M; ++m) {
      lp += arma::accu(-0.5 * LOG_2PI - arma::log(d.lam_scale) -
        0.5 * arma::square((Lam.col(m) - d.lam_loc) / d.lam_scale));
    }
    lp += arma::accu(d.psi_shape * arma::log(d.psi_scale) -
                     std::lgamma(d.psi_shape) - (d.psi_shape + 1.0) * x -
                     d.psi_scale / psi);
    lp += arma::accu(x);
  }

  if (d.scale) {
    arma::vec x = u.subvec(pos, pos + J - 1);
    pos += J;
    arma::vec sigma = arma::exp(x);
    const double df = d.sig_df;
    lp += arma::accu(std::lgamma((df + 1.0) / 2.0) - std::lgamma(df / 2.0) -
                     0.5 * std::log(df * M_PI) - std::log(d.sig_scale) -
                     (df + 1.0) / 2.0 *
                       arma::log1p(arma::square((sigma - d.sig_loc) / d.sig_scale) / df));
    lp -= J * R::pt((0.0 - d.sig_loc) / d.sig_scale, df, 0, 1);
    lp += arma::accu(x);
    A.each_col() %= sigma;
    un %= arma::square(sigma);
  }

  arma::vec nu(J, arma::fill::zeros);
  if (d.raw) {
    nu = u.subvec(pos, pos + J - 1);
    pos += J;
    lp += arma::accu(-0.5 * LOG_2PI - std::log(d.nu_scale) -
                     0.5 * arma::square((nu - d.nu_loc) / d.nu_scale));
  }

  if (!d.prior_only) lp += woodbury_loglik(d, A, un, nu);
  return std::isfinite(lp) ? lp : NEG_INF;
}

// The canonical loadings T and the log Jacobian of the map from the
// coordinates phi, which includes the Jacobian of Z = T Q.
arma::mat coords_to_T(const BridgeMap& m, const arma::vec& phi, double& log_jac) {
  const int p = m.p, M = m.M;
  arma::mat T(p, M, arma::fill::zeros);
  if (m.type == 1) {
    // Cartesian: log t_ii, then the entries below the diagonal of the anchor
    // rows, then the rows past M
    arma::uword pos = M;
    log_jac = 0.0;
    for (int i = 0; i < M; ++i) {
      T(i, i) = std::exp(phi(i));
      log_jac += (M - i) * phi(i);
      for (int j = 0; j < i; ++j) T(i, j) = phi(pos++);
    }
    for (int i = M; i < p; ++i) {
      for (int j = 0; j < M; ++j) T(i, j) = phi(pos++);
    }
    return T;
  }

  // Polar: log norms of the rows, then charts of the anchor rows (on
  // hemispheres, around the pole) and of the rows past M (on full spheres,
  // around their mean direction)
  arma::vec r = arma::exp(phi.head(p));
  log_jac = M * arma::accu(phi.head(p));
  arma::mat tau(p, M, arma::fill::zeros);
  tau(0, 0) = 1.0;
  arma::uword pos = p;
  for (int i = 1; i < M; ++i) {
    arma::vec z = phi.subvec(pos, pos + i - 1);
    pos += i;
    double zeta = arma::norm(z);
    double theta = M_PI / 2.0 * std::tanh(zeta / 2.0);
    arma::vec v = zeta > 0 ? arma::vec(z / zeta) : arma::vec(i, arma::fill::zeros);
    for (int j = 0; j < i; ++j) tau(i, j) = std::sin(theta) * v(j);
    tau(i, i) = std::cos(theta);
    log_jac += (M - i - 1) * std::log(std::cos(theta)) +
      chart_log_jac(zeta, theta, i, M_PI / 2.0);
  }
  for (int i = M; i < p; ++i) {
    arma::vec z = phi.subvec(pos, pos + M - 2);
    pos += M - 1;
    double zeta = arma::norm(z);
    double theta = M_PI * std::tanh(zeta / 2.0);
    arma::vec v = zeta > 0 ? arma::vec(z / zeta) : arma::vec(M - 1, arma::fill::zeros);
    tau.row(i) = (std::cos(theta) * m.centers[i] +
      std::sin(theta) * (m.bases[i] * v)).t();
    log_jac += chart_log_jac(zeta, theta, M - 1, M_PI);
  }
  return tau.each_col() % r;
}

}  // namespace

// [[Rcpp::export]]
SEXP befa_density_cpp(const Rcpp::List& model, const Rcpp::List& data) {
  Density* d = new Density();
  d->loadings = Rcpp::as<int>(model["loadings"]);
  d->scale = Rcpp::as<bool>(model["scale"]);
  d->raw = Rcpp::as<bool>(model["raw"]);
  d->prior_only = Rcpp::as<bool>(model["prior_only"]);
  d->J = Rcpp::as<int>(model["J"]);
  d->M = Rcpp::as<int>(model["M"]);
  d->h2_a = Rcpp::as<double>(model["h2_a"]);
  d->h2_b = Rcpp::as<double>(model["h2_b"]);
  d->xi = Rcpp::as<double>(model["xi"]);
  d->log_C = Rcpp::as<double>(model["log_C"]);
  d->lam_a = Rcpp::as<double>(model["lam_a"]);
  d->lam_b = Rcpp::as<double>(model["lam_b"]);
  d->lam_loc = Rcpp::as<arma::vec>(model["lam_loc"]);
  d->lam_scale = Rcpp::as<arma::vec>(model["lam_scale"]);
  d->psi_shape = Rcpp::as<double>(model["psi_shape"]);
  d->psi_scale = Rcpp::as<arma::vec>(model["psi_scale"]);
  d->sig_df = Rcpp::as<double>(model["sig_df"]);
  d->sig_loc = Rcpp::as<double>(model["sig_loc"]);
  d->sig_scale = Rcpp::as<double>(model["sig_scale"]);
  d->nu_loc = Rcpp::as<double>(model["nu_loc"]);
  d->nu_scale = Rcpp::as<double>(model["nu_scale"]);
  d->obs_mat = Rcpp::as<arma::mat>(data["obs_mat"]);
  d->root_pad = Rcpp::as<arma::mat>(data["root_pad"]);
  d->n_obs = Rcpp::as<arma::vec>(data["n_obs"]);
  d->r2_item = Rcpp::as<arma::vec>(data["r2_item"]);
  d->ybar_pad = Rcpp::as<arma::mat>(data["ybar_pad"]);
  d->n_p = Rcpp::as<arma::vec>(data["n_p"]);
  d->K_p = Rcpp::as<arma::uvec>(data["K_p"]);
  d->ll_const = Rcpp::as<double>(data["ll_const"]);
  return Rcpp::XPtr<Density>(d, true);
}

// [[Rcpp::export]]
double befa_log_density_cpp(SEXP density, const arma::vec& u) {
  Rcpp::XPtr<Density> d(density);
  return log_density(*d, u);
}

// [[Rcpp::export]]
SEXP bridge_map_cpp(const Rcpp::List& map) {
  BridgeMap* m = new BridgeMap();
  m->type = Rcpp::as<int>(map["type"]);
  m->p = Rcpp::as<int>(map["p"]);
  m->M = Rcpp::as<int>(map["M"]);
  m->n_u = Rcpp::as<int>(map["n_u"]);
  m->inv_ord = Rcpp::as<arma::uvec>(map["inv_ord"]) - 1;
  m->load_idx = Rcpp::as<arma::uvec>(map["load_idx"]) - 1;
  m->other_idx = Rcpp::as<arma::uvec>(map["other_idx"]) - 1;
  m->log_vol = Rcpp::as<double>(map["log_vol"]);
  Rcpp::List centers = map["centers"], bases = map["bases"];
  m->centers.resize(m->p);
  m->bases.resize(m->p);
  for (int i = m->M; i < m->p && m->type == 0; ++i) {
    m->centers[i] = Rcpp::as<arma::vec>(centers[i]);
    m->bases[i] = Rcpp::as<arma::mat>(bases[i]);
  }
  return Rcpp::XPtr<BridgeMap>(m, true);
}

// [[Rcpp::export]]
arma::mat bridge_coords_to_T_cpp(SEXP map, const arma::vec& phi) {
  Rcpp::XPtr<BridgeMap> m(map);
  double log_jac = 0.0;
  return coords_to_T(*m, phi, log_jac);
}

// [[Rcpp::export]]
double bridge_log_posterior_cpp(SEXP density, SEXP map, const arma::vec& x) {
  Rcpp::XPtr<Density> d(density);
  Rcpp::XPtr<BridgeMap> m(map);
  const arma::uword n_other = m->other_idx.n_elem;
  arma::vec u(m->n_u);
  if (n_other > 0) u.elem(m->other_idx) = x.head(n_other);
  double log_jac = 0.0;
  arma::mat T = coords_to_T(*m, x.tail(x.n_elem - n_other), log_jac);
  arma::mat Z = T.rows(m->inv_ord);
  u.elem(m->load_idx) = arma::vectorise(Z);
  double lp = log_density(*d, u) + log_jac + m->log_vol;
  return std::isfinite(lp) ? lp : NEG_INF;
}
