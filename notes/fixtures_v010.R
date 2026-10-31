# fixtures_v010.R
# Build the golden test fixtures from BayesianEFA v0.1.0 (commit 24c15ac)
# Author: Ricardo Rey-Sáez
# Last modified: 31-10-2026

# Run it from the package root with `Rscript notes/fixtures_v010.R`. It writes
# tests/testthat/fixtures/*_v010.rds, each with fixed inputs and the output
# v0.1.0 returns for them. The v0.1.0 sources are read from git, so the script
# keeps working after the rewrite. Only synthetic data are stored.

base_commit <- "24c15ac"
src_dir <- file.path(tempdir(), "befa_v010")
out_dir <- file.path("tests", "testthat", "fixtures")
dir.create(src_dir, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (pkg in c("Rcpp", "RcppArmadillo", "posterior", "mvnfast")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' is required.", call. = FALSE)
}

# Every fixture records where and when it was made.
fixture_meta <- list(
  commit = base_commit,
  created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  r_version = R.version.string
)

# v0.1.0 sources

# This function writes a file of the base commit, byte by byte, to the
# temporary folder and returns its path.
git_show <- function(path) {
  dest <- file.path(src_dir, basename(path))
  status <- system2("git", c("show", paste0(base_commit, ":", path)), stdout = dest)
  if (!identical(status, 0L)) stop("git show failed for ", path, call. = FALSE)
  dest
}

Rcpp::sourceCpp(git_show("src/rsp_efficient.cpp"))
Rcpp::sourceCpp(git_show("src/factor_scores_cpp.cpp"))
Rcpp::sourceCpp(git_show("src/posterior_metrics_cpp.cpp"))
source(git_show("R/rsp_align.R"))
source(git_show("R/utils_FIML.R"))

# This function saves one fixture as list(inputs, output, meta).
save_fixture <- function(name, inputs, output) {
  file <- file.path(out_dir, paste0(name, "_v010.rds"))
  saveRDS(list(inputs = inputs, output = output, meta = fixture_meta), file)
  cat(sprintf("  %-40s %6.1f KB\n", basename(file), file.size(file) / 1024))
}

# 1. rsp_align()

# This function returns draws of a simple-structure loading matrix, each one
# rotated by a random orthogonal matrix, which includes reflections and
# permutations.
make_rsp_draws <- function(J, M, n_chains, n_iter, seed) {
  set.seed(seed)
  Lambda_true <- matrix(0.1, nrow = J, ncol = M)
  for (j in seq_len(J)) Lambda_true[j, (j - 1) %% M + 1] <- 0.7
  draws <- matrix(NA_real_, nrow = n_chains * n_iter, ncol = J * M)
  for (s in seq_len(nrow(draws))) {
    L_s <- Lambda_true + matrix(rnorm(J * M, 0, 0.05), nrow = J, ncol = M)
    Q_s <- qr.Q(qr(matrix(rnorm(M * M), nrow = M, ncol = M)))
    draws[s, ] <- c(L_s %*% Q_s)
  }
  draws
}

cat("rsp_align():\n")
rsp_cases <- list(
  list(name = "rsp_align_J9_M3", J = 9, M = 3, n_chains = 4, n_iter = 60, seed = 101),
  list(name = "rsp_align_J6_M2", J = 6, M = 2, n_chains = 2, n_iter = 100, seed = 102),
  list(name = "rsp_align_J5_M1", J = 5, M = 1, n_chains = 2, n_iter = 100, seed = 103)
)
for (cs in rsp_cases) {
  draws <- make_rsp_draws(cs$J, cs$M, cs$n_chains, cs$n_iter, cs$seed)
  for (fmt in c("column_major", "row_major")) {
    input <- draws
    if (fmt == "row_major") {
      # The same draws, flattened item by item.
      input <- draws[, as.vector(t(matrix(seq_len(cs$J * cs$M), nrow = cs$J, ncol = cs$M)))]
    }
    res <- rsp_align(
      lambda_draws = input, n_items = cs$J, n_factors = cs$M,
      n_chains = cs$n_chains, format = fmt
    )
    # The summary depends on the posterior version, not on v0.1.0.
    res$summary <- NULL
    save_fixture(
      paste0(cs$name, "_", fmt),
      inputs = list(lambda_draws = input, n_items = cs$J, n_factors = cs$M,
                    n_chains = cs$n_chains, format = fmt),
      output = res
    )
  }
}

# 2. compute_factor_scores_cpp()

# This function returns correlation-model draws laid out as v0.1.0 reads them:
# Lambda (column-major) followed by Rho.
make_cor_draws <- function(J, M, n_iter, n_chains, seed) {
  set.seed(seed)
  Lambda_true <- matrix(0.15, nrow = J, ncol = M)
  for (j in seq_len(J)) Lambda_true[j, (j - 1) %% M + 1] <- 0.65
  n_vars <- J * M + J * J
  cube <- array(NA_real_, dim = c(n_iter, n_chains, n_vars))
  for (ch in seq_len(n_chains)) {
    for (s in seq_len(n_iter)) {
      L_s <- Lambda_true + matrix(rnorm(J * M, 0, 0.03), nrow = J, ncol = M)
      Rho_s <- tcrossprod(L_s)
      diag(Rho_s) <- 1
      cube[s, ch, ] <- c(c(L_s), c(Rho_s))
    }
  }
  list(cube = cube, Lambda_true = Lambda_true)
}

cat("compute_factor_scores_cpp():\n")
J <- 6; M <- 2; N <- 40
cor_draws <- make_cor_draws(J, M, n_iter = 20, n_chains = 2, seed = 201)
Rho_true <- tcrossprod(cor_draws$Lambda_true)
diag(Rho_true) <- 1
set.seed(202)
Y <- scale(mvnfast::rmvn(N, mu = rep(0, J), sigma = Rho_true))
Y_na <- Y
Y_na[sample(length(Y_na), 25)] <- NA

for (cs in list(list(name = "factor_scores_cor_complete", data = unclass(Y)),
                list(name = "factor_scores_cor_missing", data = unclass(Y_na)))) {
  data <- cs$data
  attributes(data) <- list(dim = dim(data))
  obs_patterns_0 <- lapply(seq_len(N), function(i) as.integer(which(!is.na(data[i, ])) - 1L))
  inputs <- list(
    post_draws = cor_draws$cube, data = data,
    lambda_idx = as.integer(0:(J * M - 1)),
    sig_idx = as.integer(J * M + 0:(J * J - 1)),
    nu_idx = integer(0), J = J, M = M, N = N,
    has_missing = anyNA(data), obs_patterns = obs_patterns_0,
    seed = 203
  )
  set.seed(inputs$seed)
  eta <- compute_factor_scores_cpp(
    post_draws = inputs$post_draws, data = inputs$data,
    lambda_idx = inputs$lambda_idx, sig_idx = inputs$sig_idx, nu_idx = inputs$nu_idx,
    J = J, M = M, N = N, has_missing = inputs$has_missing,
    obs_patterns = inputs$obs_patterns
  )
  save_fixture(cs$name, inputs = inputs, output = eta)
}

# 3. compute_posterior_metrics_cpp()

cat("compute_posterior_metrics_cpp():\n")

# This function repeats the steps of compute_posterior_metrics() in v0.1.0 for
# the given draws, and returns the inputs it passed to C++ and the output.
run_metrics <- function(cube, data, model_type, seed) {
  N <- nrow(data)
  J <- ncol(data)
  has_missing <- anyNA(data)
  comp_idx <- which(stats::complete.cases(data))
  miss_idx <- which(!stats::complete.cases(data))

  if (has_missing) {
    h1 <- estimate_h1_moments_em(data)
    M_sat <- h1$Mu
    S_sat <- h1$Sigma
    if (model_type == 3) {
      D_inv <- diag(1 / sqrt(diag(S_sat)))
      S_sat <- D_inv %*% S_sat %*% D_inv
      M_sat <- rep(0, J)
    }
    sat_fiml <- list(comp_idx = comp_idx, miss_idx = miss_idx,
                     Y_comp = data[comp_idx, , drop = FALSE],
                     Y_miss = data[miss_idx, , drop = FALSE])
    ll_saturated <- sum(compute_fiml_loglik(sat_fiml, mu = M_sat, Sigma = S_sat))
  } else {
    S_sat <- if (model_type == 3) stats::cor(data) else stats::cov(data) * (N - 1) / N
    M_sat <- if (model_type == 1) colMeans(data) else rep(0, J)
    ll_saturated <- sum(mvnfast::dmvn(data, mu = M_sat, sigma = S_sat, log = TRUE))
  }

  n_sig <- J * J
  inputs <- list(
    post_draws = cube, data = data,
    sig_idx = as.integer(0:(n_sig - 1)),
    nu_idx = if (model_type == 1) as.integer(n_sig + 0:(J - 1)) else integer(0),
    M_sat = M_sat, S_sat = S_sat, ll_saturated = ll_saturated,
    model_type = as.integer(model_type), has_missing = has_missing,
    Y_comp = if (length(comp_idx)) data[comp_idx, , drop = FALSE] else matrix(0, 0, J),
    Y_miss = if (length(miss_idx)) data[miss_idx, , drop = FALSE] else matrix(0, 0, J),
    comp_idx_0 = as.integer(comp_idx - 1L), miss_idx_0 = as.integer(miss_idx - 1L),
    obs_patterns_miss = lapply(miss_idx, function(i) as.integer(which(!is.na(data[i, ])) - 1L)),
    seed = seed
  )
  set.seed(seed)
  output <- do.call(compute_posterior_metrics_cpp, inputs[setdiff(names(inputs), "seed")])
  list(inputs = inputs, output = output)
}

# This function returns raw-model draws laid out as v0.1.0 reads them: Sigma
# (column-major) followed by nu.
make_raw_draws <- function(Lambda, psi, mu, n_iter, n_chains, seed) {
  set.seed(seed)
  J <- nrow(Lambda)
  cube <- array(NA_real_, dim = c(n_iter, n_chains, J * J + J))
  for (ch in seq_len(n_chains)) {
    for (s in seq_len(n_iter)) {
      L_s <- Lambda + matrix(rnorm(length(Lambda), 0, 0.03), nrow = J)
      Sigma_s <- tcrossprod(L_s) + diag(psi * exp(rnorm(J, 0, 0.05)))
      cube[s, ch, ] <- c(c(Sigma_s), mu + rnorm(J, 0, 0.05))
    }
  }
  cube
}

set.seed(301)
Lambda_raw <- cor_draws$Lambda_true * 2
psi_raw <- rep(1.5, J)
mu_raw <- seq(1, 3, length.out = J)
Y_raw <- mvnfast::rmvn(N, mu = mu_raw, sigma = tcrossprod(Lambda_raw) + diag(psi_raw))
Y_raw_na <- Y_raw
Y_raw_na[sample(length(Y_raw_na), 25)] <- NA
raw_cube <- make_raw_draws(Lambda_raw, psi_raw, mu_raw, n_iter = 20, n_chains = 2, seed = 302)

# The correlation model reads only Rho, so the Lambda block is dropped.
rho_cube <- cor_draws$cube[, , J * M + seq_len(J * J), drop = FALSE]

metric_cases <- list(
  list(name = "posterior_metrics_cor_complete", cube = rho_cube, data = unclass(Y), type = 3),
  list(name = "posterior_metrics_raw_complete", cube = raw_cube, data = Y_raw, type = 1),
  list(name = "posterior_metrics_raw_missing", cube = raw_cube, data = Y_raw_na, type = 1)
)
for (cs in metric_cases) {
  data <- cs$data
  attributes(data) <- list(dim = dim(data))
  res <- run_metrics(cs$cube, data, cs$type, seed = 303)
  save_fixture(cs$name, inputs = res$inputs, output = res$output)
}

# 4. FIML log-likelihood and EM saturated moments (R code in v0.1.0)

cat("compute_fiml_loglik() and estimate_h1_moments_em():\n")
data <- Y_raw_na
comp_idx <- which(stats::complete.cases(data))
miss_idx <- which(!stats::complete.cases(data))
data_fiml <- list(comp_idx = comp_idx, miss_idx = miss_idx,
                  Y_comp = data[comp_idx, , drop = FALSE],
                  Y_miss = data[miss_idx, , drop = FALSE])
sigma_null <- apply(data, 2, sd, na.rm = TRUE)
mu_null <- colMeans(data, na.rm = TRUE)

save_fixture(
  "fiml_loglik_null_diag",
  inputs = list(data_fiml = data_fiml, mu = mu_null, Sigma = diag(sigma_null^2)),
  output = compute_fiml_loglik(data_fiml, mu = mu_null, Sigma = diag(sigma_null^2))
)
save_fixture(
  "em_h1_moments",
  inputs = list(Y = data),
  output = estimate_h1_moments_em(data)
)

cat("Done.\n")
