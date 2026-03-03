#' Bayesian Factor Scores
#'
#' Computes posterior factor scores from a fitted BEFA model by evaluating the
#' conditional posterior distribution of the latent factors given the observed data
#' and the estimated model parameters.
#'
#' @param object A `befa` object returned by [befa()].
#' @param post_summaries Logical. If `TRUE` (default), includes posterior summaries
#'   computed via [posterior::summarise_draws()]. Set to `FALSE` to skip summaries.
#'
#' @details
#' **Factor Score Estimation**
#'
#' Consistent with the data generation process of the factor model, the observation
#' for the \eqn{i}-th individual is defined as:
#' \deqn{\mathbf{y}_i = \boldsymbol{\nu} + \mathbf{\Lambda} \cdot \boldsymbol{\eta}_i + \boldsymbol{\epsilon}_i}
#'
#' where \eqn{\boldsymbol{\eta}_i \sim \mathcal{N}(\mathbf{0}, \mathbf{I}_M)} and
#' \eqn{\boldsymbol{\epsilon}_i \sim \mathcal{N}(\mathbf{0}, \mathbf{\Psi})}.
#'
#' Because the posterior distributions of the model parameters are available from
#' the MCMC draws, the factor scores can be derived analytically for each iteration
#' (see, e.g., Levy & Mislevy, 2016, pp.199--200). The conditional posterior
#' distribution of the latent factors given the observed data is:
#' \deqn{\boldsymbol{\eta}_i \mid \mathbf{y}_i, \mathbf{\Lambda}, \mathbf{\Psi}, \boldsymbol{\nu} \sim \mathcal{N}(\boldsymbol{\mu}_{\eta \mid y}, \mathbf{\Sigma}_{\eta \mid y})}
#'
#' with the conditional covariance matrix and mean vector defined as:
#' \deqn{\mathbf{\Sigma}_{\eta \mid y} = (\mathbf{I}_M + \mathbf{\Lambda}^\top \cdot \mathbf{\Psi}^{-1} \cdot \mathbf{\Lambda})^{-1}}
#' \deqn{\boldsymbol{\mu}_{\eta \mid y} = \mathbf{\Sigma}_{\eta \mid y} \cdot \mathbf{\Lambda}^\top \cdot \mathbf{\Psi}^{-1} \cdot (\mathbf{y}_i - \boldsymbol{\nu})}
#'
#' **Handling of Missing Data**
#'
#' For incomplete observations, the estimation naturally adapts by conditioning
#' only on the observed subset of variables for person \eqn{i}. Specifically, the
#' matrices \eqn{\mathbf{\Lambda}} and \eqn{\mathbf{\Psi}}, as well as the
#' data vector \eqn{(\mathbf{y}_i - \boldsymbol{\nu})}, are subsetted to the
#' indices of the observed variables, denoted as \eqn{o_i}.
#' @return A list of class `befa_factor_scores` containing:
#'   * `eta_draws`: A `draws_array` (iterations \eqn{\times} chains \eqn{\times} P)
#'     containing the posterior draws of the factor scores.
#'   * `summary`: (If `post_summaries = TRUE`) A `draws_summary` data frame computed
#'     via [posterior::summarise_draws()].
#'
#' @references
#' Levy, R., & Mislevy, R. J. (2016). *Bayesian psychometric modeling*.
#' Chapman & Hall/CRC.
#'
#' @examples
#' \dontrun{
#' # ------------------------------------------------------- #
#' #    1. Factor scores after fitting Bayesian EFA model    #
#' # ------------------------------------------------------- #
#'
#' # Fit Bayesian EFA model to the famous Grant-White School Data (Holzinger y Swineford , 1939)
#' befa_fit <- befa(
#'   data = HS_data,
#'   n_factors = 3,
#'   factor_scores = FALSE,
#'   compute_fit_indices = FALSE,
#'   compute_reliability = FALSE,
#'   backend = "rstan",
#'   seed = 17,
#'   chains = 4,
#'   parallel_chains = 4
#' )
#'
#' # Compute factor scores from a fitted befa object
#' Bayes_fsc <- befa_factor_scores(befa_fit, post_summaries = TRUE)
#'
#' # Store Posterior Means in a 300 x 2 Matrix
#' Bayes_fsc_pmean <- matrix(Bayes_fsc$summary$mean, nrow = nrow(HS_data), ncol = 3)
#'
#' # Estimate frequentist EFA and store factor scores
#' efa_fit <- factanal(x = HS_data, factors = 3, scores = "regression")
#'
#' # Equivalence between factor scores (we order the columns)
#' cor(efa_fit$scores[, 2], Bayes_fsc_pmean[, 1])
#' cor(efa_fit$scores[, 1], Bayes_fsc_pmean[, 2])
#' cor(efa_fit$scores[, 3], Bayes_fsc_pmean[, 3])
#'
#' # Bivariate plot to see linear equivalence
#' par(mfrow = c(1, 3))
#' plot(Bayes_fsc_pmean[, 1], efa_fit$scores[, 2],
#'      main = "Factor 1",
#'      xlab = "Bayesian factor scores",
#'      ylab = "Frequentist factor scores"
#' )
#' abline(a = 0, b = 1, col = "steelblue2", lwd = 2)
#' plot(Bayes_fsc_pmean[, 2], efa_fit$scores[, 1],
#'      main = "Factor 2",
#'      xlab = "Bayesian factor scores",
#'      ylab = "Frequentist factor scores"
#' )
#' abline(a = 0, b = 1, col = "salmon2", lwd = 2)
#' plot(Bayes_fsc_pmean[, 3], efa_fit$scores[, 3],
#'      main = "Factor 3",
#'      xlab = "Bayesian factor scores",
#'      ylab = "Frequentist factor scores"
#' )
#' abline(a = 0, b = 1, col = "firebrick4", lwd = 2)
#' par(mfrow = c(1, 1))
#'
#' # -----------------------------------------------------------------------------
#'
#' # ------------------------------------------------------ #
#' #    2. Factor scores when fitting Bayesian EFA model    #
#' # ------------------------------------------------------ #
#' # Fit Bayesian model with Factor Scores = TRUE
#' befa_fit <- befa(
#'   data = HS_data,
#'   n_factors = 3,
#'   factor_scores = TRUE,
#'   compute_fit_indices = FALSE,
#'   compute_reliability = FALSE,
#'   backend = "rstan",
#'   seed = 17,
#'   chains = 4,
#'   parallel_chains = 4
#' )
#'
#' # Factor scores (eta) are stored inside the befa object
#' extract_posterior_draws(befa_fit, pars = "eta")
#'
#' # # A draws_matrix: 1000 iterations, 4 chains, and 903 variables
#' #     variable
#' # draw eta[1,1] eta[2,1] eta[3,1] eta[4,1] eta[5,1] eta[6,1] eta[7,1] eta[8,1]
#' #   1    -1.215     0.55   -0.437    0.158    0.716     0.77    -1.23     0.78
#' #   2     0.098    -0.43   -0.645    1.588   -0.511     1.21    -1.15    -0.66
#' #   3    -1.324     0.36    0.116    0.363   -1.512     0.15    -1.73     0.14
#' #   4     0.578     1.79   -0.108   -0.024    0.054     1.44    -0.61    -1.55
#' #   5    -1.244     1.06   -0.393    1.662   -0.938     0.28    -0.52    -0.55
#' #   6    -1.569     1.23   -0.059    0.988   -1.944     0.97    -1.29     0.53
#' #   7    -1.099     0.35   -0.372    1.401    0.407     0.40    -0.50     0.55
#' #   8     0.043     1.40   -0.198    0.997   -0.279     1.59    -1.15    -0.16
#' #   9     0.078     0.60    0.129    0.833   -1.312     0.90    -0.41     0.31
#' #   10   -0.084     0.72   -0.665    0.573   -1.653     0.99    -0.95    -0.72
#' # # ... with 3990 more draws, and 895 more variables
#'
#' posterior_summaries(befa_fit, pars = "eta")
#'
#' # # A tibble: 903 × 10
#' #    variable    mean median    sd   mad     q5    q95  rhat ess_bulk ess_tail
#' #    <chr>      <dbl>  <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>    <dbl>    <dbl>
#' #  1 eta[1,1]  -0.719 -0.734 0.627 0.645 -1.73   0.317 1.00     4025.    3803.
#' #  2 eta[2,1]   0.595  0.591 0.624 0.611 -0.430  1.61  1.00     4003.    3849.
#' #  3 eta[3,1]  -0.161 -0.165 0.601 0.588 -1.16   0.836 1.00     3732.    4035.
#' #  4 eta[4,1]   0.717  0.706 0.600 0.594 -0.264  1.72  1.00     3851.    3791.
#' #  5 eta[5,1]  -0.507 -0.507 0.596 0.586 -1.47   0.497 1.000    3925.    3645.
#' #  6 eta[6,1]   0.469  0.473 0.611 0.624 -0.528  1.48  1.000    3556.    3836.
#' #  7 eta[7,1]  -1.40  -1.40  0.611 0.622 -2.40  -0.422 1.00     3698.    3830.
#' #  8 eta[8,1]  -0.134 -0.127 0.596 0.580 -1.11   0.846 1.00     4070.    3815.
#' #  9 eta[9,1]  -0.298 -0.298 0.595 0.600 -1.28   0.672 1.000    4000.    3938.
#' # 10 eta[10,1] -1.35  -1.35  0.593 0.605 -2.29  -0.358 1.00     3965.    3718.
#' # # 893 more rows
#' # # Use `print(n = ...)` to see more rows
#'
#' # You can store posterior means easely because Stan always use column-major order
#' Bayes_fsc <- matrix(
#'   posterior_summaries(befa_fit, pars = "eta")$mean,
#'   nrow = nrow(HS_data),
#'   ncol = 3
#' )
#'
#' # -----------------------------------------------------------------------------
#' }
#'
#' @export
befa_factor_scores <- function(object, post_summaries = TRUE) {
  # ────────────────────── #
  #    Input validation    #
  # ────────────────────── #

  if (!inherits(object, "befa")) {
    stop("Input 'object' must be of class 'befa'.", call. = FALSE)
  }

  # Raw data is required for factor scores
  if (!is.null(object$stan_data$Y_original)) {
    data <- as.matrix(object$stan_data$Y_original)
  } else if (!is.null(object$stan_data$Y) && nrow(object$stan_data$Y) > 0) {
    data <- as.matrix(object$stan_data$Y)
  } else {
    stop("Factor scores require raw data. Cannot compute from summary statistics.",
      call. = FALSE
    )
  }

  # ─────────────────────────── #
  #    Filter and scale data    #
  # ─────────────────────────── #

  # Drop observations the model did not use (e.g., listwise deletion for cor/cov)
  N_model <- object$stan_data$N_complete + object$stan_data$N_incomplete
  if (nrow(data) > N_model) {
    data <- data[complete.cases(data), , drop = FALSE]
  }

  # Scale data to match model parameterization:
  model_type <- object$model_type
  if (model_type == "cor") {
    col_means <- colMeans(data, na.rm = TRUE)
    col_sds <- apply(data, 2, sd, na.rm = TRUE)
    data <- scale(data, center = col_means, scale = col_sds)
  } else if (model_type == "cov") {
    col_means <- colMeans(data, na.rm = TRUE)
    data <- scale(data, center = col_means, scale = FALSE)
  }

  N <- nrow(data)
  J <- ncol(data)
  M <- object$n_factors

  # ───────────────────────────── #
  #    Extract posterior draws    #
  # ───────────────────────────── #

  # Use as_draws_array to preserve chain structure: (iterations x chains x variables)
  post_draws <- posterior::as_draws_array(object$stanfit)
  n_iter <- dim(post_draws)[1]
  n_chains <- dim(post_draws)[2]
  var_names <- dimnames(post_draws)[[3]]

  # Parameter indices
  lambda_idx <- grep("^Lambda\\[", var_names)

  if (model_type == "raw") {
    sig_idx <- grep("^Sigma\\[", var_names)
    nu_idx <- grep("^nu\\[", var_names)
  } else if (model_type == "cov") {
    sig_idx <- grep("^Sigma\\[", var_names)
    nu_idx <- NULL
  } else { # cor
    sig_idx <- grep("^Rho\\[", var_names)
    nu_idx <- NULL
  }

  # ────────────────────────────────────── #
  #    Precompute missing data patterns    #
  # ────────────────────────────────────── #

  obs_patterns <- lapply(seq_len(N), function(i) which(!is.na(data[i, ])))
  has_missing <- anyNA(data)

  # ──────────────────────────────────── #
  #    Prepare inputs for C++ backend    #
  # ──────────────────────────────────── #

  # Convert to 0-based indices for C++
  lambda_idx_0 <- as.integer(lambda_idx - 1L)
  sig_idx_0 <- as.integer(sig_idx - 1L)
  nu_idx_0 <- if (!is.null(nu_idx)) as.integer(nu_idx - 1L) else integer(0)

  # Observation patterns as 0-based integer vectors
  obs_patterns_0 <- lapply(obs_patterns, function(x) as.integer(x - 1L))

  # Convert draws to a 3D array (n_iter x n_chains x n_vars) for C++
  draws_cube <- as.array(post_draws)

  # ────────────────────────────────── #
  #    Compute factor scores in C++    #
  # ────────────────────────────────── #

  eta_cube <- compute_factor_scores_cpp(
    post_draws = draws_cube,
    data = data,
    lambda_idx = lambda_idx_0,
    sig_idx = sig_idx_0,
    nu_idx = nu_idx_0,
    J = J, M = M, N = N,
    has_missing = has_missing,
    obs_patterns = obs_patterns_0
  )

  # ──────────────────────────────── #
  #    Wrap result as draws_array    #
  # ──────────────────────────────── #

  # Build Stan-style parameter names in column-major order:
  # eta[1,1], eta[2,1], ..., eta[N,1], eta[1,2], eta[2,2], ..., eta[N,M]
  eta_names <- paste0(
    "eta[",
    rep(seq_len(N), times = M), ",",
    rep(seq_len(M), each = N), "]"
  )

  dimnames(eta_cube) <- list(
    iteration = seq_len(n_iter),
    chain     = seq_len(n_chains),
    variable  = eta_names
  )

  eta_array <- eta_cube

  # ────────────────────────────────────── #
  #    Convert to posterior draws_array    #
  # ────────────────────────────────────── #

  # eta_array is already (iterations x chains x P) with proper dimnames
  eta_draws <- posterior::as_draws_array(eta_array)

  # ───────────────────────── #
  #    Build return object    #
  # ───────────────────────── #

  result <- list(
    eta_draws = eta_draws
  )

  # Compute posterior summaries if requested
  if (post_summaries) {
    result$summary <- posterior::summarise_draws(eta_draws)
  }

  class(result) <- "befa_factor_scores"
  return(result)
}


#' Inject Factor Scores into stanfit Object (Internal)
#'
#' Takes computed factor scores and injects them as additional parameters
#' into the stanfit object. N and M are derived from the parameter names.
#'
#' @param stanfit A stanfit object.
#' @param eta_draws A `draws_array` (iterations x chains x P) of factor scores.
#' @return Updated stanfit object with eta parameters.
#' @keywords internal
#' @noRd
inject_factor_scores <- function(stanfit, eta_draws) {
  sim <- stanfit@sim
  n_chains <- length(sim$samples)
  n_kept <- sim$n_save[1] - sim$warmup2[1]

  # Parameter names from the 3rd dimension of eta_draws
  eta_names <- dimnames(eta_draws)[[3]]

  # Deduce N and M from the parameter names (assuming eta[i,m] format)
  # Extract indices from "eta[i,m]"
  indices <- do.call(rbind, lapply(strsplit(gsub("eta\\[|\\]", "", eta_names), ","), as.numeric))
  N <- max(indices[, 1])
  M <- max(indices[, 2])

  # Inject into each chain (following the RSP injection pattern)
  for (ch in seq_len(n_chains)) {
    # Get the total vector length from an existing parameter (includes warmup slots)
    existing_vec <- sim$samples[[ch]][[1]]
    total_len <- length(existing_vec)
    insertion_start <- total_len - n_kept + 1

    for (p in seq_along(eta_names)) {
      # Build a full-length vector: warmup zeros + post-warmup values
      new_vec <- rep(0, total_len)
      new_vec[insertion_start:total_len] <- as.numeric(eta_draws[, ch, p])
      sim$samples[[ch]][[eta_names[p]]] <- new_vec
    }
  }

  # Update parameter metadata
  sim$pars_oi <- c(sim$pars_oi, "eta")
  sim$dims_oi[["eta"]] <- c(N, M)
  sim$fnames_oi <- c(sim$fnames_oi, eta_names)
  sim$n_flatnames <- length(sim$fnames_oi)

  stanfit@sim <- sim
  return(stanfit)
}
