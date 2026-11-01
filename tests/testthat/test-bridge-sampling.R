# Bridge sampling on the canonical representative of the rotations

# This function draws S independent values of the Stan parameters of a
# structure from their prior (two "chains" of S / 2), and returns them with
# the posterior-mean Gram matrix of the standardized loadings.
.prior_draws <- function(spec, sd, S = 4000L, seed = 1) {
  set.seed(seed)
  J <- sd$J
  M <- if (spec$loadings == "uni") 1L else sd$M
  pars <- list()
  if (spec$loadings == "uv") {
    pars$h2 <- matrix(stats::rbeta(S * J, sd$pr_h2[1], sd$pr_h2[2]), S)
    dir <- array(stats::rnorm(S * J * M), c(S, J, M))
    norm <- sqrt(apply(dir^2, c(1, 2), sum))
    r <- matrix(sqrt(stats::rchisq(S * J, sd$pr_xi + M)), S)
    pars$Z <- matrix(sweep(dir, c(1, 2), r / norm, `*`), S)
    L_std <- matrix(sweep(dir, c(1, 2), sqrt(pars$h2) / norm, `*`), S)
  } else if (spec$loadings == "uni") {
    pars$lambda_std <- matrix(2 * stats::rbeta(S * J, sd$pr_lambda[1], sd$pr_lambda[2]) - 1, S)
    L_std <- pars$lambda_std
  } else {
    pars$Lambda <- matrix(stats::rnorm(S * J * M, rep(rep(sd$pr_lambda_loc, M), each = S),
                                       rep(rep(sd$pr_lambda_scale, M), each = S)), S)
    pars$psi <- matrix(1 / stats::rgamma(S * J, sd$pr_psi_shape,
                                         rate = rep(sd$pr_psi_scale, each = S)), S)
    variance <- .item_sums(pars$Lambda^2, J, M) + pars$psi
    L_std <- pars$Lambda / .repeat_items(sqrt(variance), M)
  }
  if (.spec_has_scale(spec)) {
    pars$sigma <- matrix(abs(stats::rt(S * J, sd$pr_sigma[1])) * sd$pr_sigma[3], S)
  }
  if (spec$model_type == "raw") {
    pars$nu <- matrix(stats::rnorm(S * J, sd$pr_nu[1], sd$pr_nu[2]), S)
  }
  first <- rep(seq_len(S / 2) <= S / 4, 2)
  G <- Reduce(`+`, lapply(which(first), function(s) tcrossprod(matrix(L_std[s, ], J, M))))
  list(pars = pars, G = G / sum(first))
}

# This function estimates the log marginal likelihood of the prior alone,
# which is exactly 0.
.prior_logml <- function(spec, J = 4L, M = 2L, anchors = NULL) {
  prep <- .prepare_befa_data(.simulated_data(N = 40, J = 5)[, seq_len(J)],
                             spec$model_type)
  prior <- prepare_befa_priors(list(), spec$model_type,
                               if (spec$loadings == "normal") "normal" else "unit_vector", M)
  sd <- .stan_data(spec, prep, M, prior)
  draws <- .prior_draws(spec, sd)
  setup <- .bridge_setup(spec, sd, draws$pars, n_chains = 2L, G = draws$G,
                         anchors = anchors, prior_only = TRUE)
  res <- bridgesampling::bridge_sampler(
    samples = setup$samples, log_posterior = .bridge_log_posterior,
    data = setup$data, lb = setup$lb, ub = setup$ub, method = "warp3",
    silent = TRUE)
  list(logml = res$logml, setup = setup)
}

test_that("the volume of the orthogonal group is that of Muirhead (1982)", {
  expect_equal(.log_vol_orthogonal(1), log(2))
  expect_equal(.log_vol_orthogonal(2), log(4 * pi))
  expect_equal(.log_vol_orthogonal(3), log(16 * pi^2))
})

test_that("the anchors are strong items of distinct factors", {
  # Item 5 has the largest communality; item 1, the largest one left once
  # item 5 is explained.
  L <- matrix(c(0.8, 0.7, 0.75, 0.1, 0.1, 0.1, 0.1, 0.1, 0.6, 0.9), 5, 2)
  piv <- .select_anchors(tcrossprod(L), 2L)
  expect_identical(piv$anchors, c(5L, 1L))
  expect_equal(piv$pivots[1], sum(L[5, ]^2))
  expect_identical(.select_anchors(tcrossprod(L), 2L, fixed = c(2L, 4L))$anchors, c(2L, 4L))
  expect_error(.select_anchors(tcrossprod(L), 2L, fixed = c(2L, 2L)), "2 distinct")
})

test_that("the canonical T is lower trapezoidal and the same for every rotation", {
  set.seed(3)
  Z <- matrix(stats::rnorm(15), 5, 3)
  ord <- c(2L, 5L, 1L, 3L, 4L)
  T <- .rotate_to_T(Z, ord, 3L)
  anchor_block <- T[1:3, ]
  expect_true(all(abs(anchor_block[upper.tri(anchor_block)]) < 1e-12))
  expect_true(all(diag(anchor_block) > 0))
  expect_equal(tcrossprod(T), tcrossprod(Z[ord, ]))
  Q <- qr.Q(qr(matrix(stats::rnorm(9), 3)))
  expect_equal(.rotate_to_T(Z %*% Q, ord, 3L), T)
  # One factor: the sign of the anchor
  z <- matrix(c(-0.4, 0.3, 0.9), 3)
  expect_equal(.rotate_to_T(z, c(1L, 2L, 3L), 1L), -z)
})

test_that("bridge sampling recovers the evidence 1 of the prior", {
  skip_if_not_installed("bridgesampling")
  skip_if_not_installed("coda")
  # Unit-vector prior (polar coordinates), with scale and means
  uv <- .prior_logml(.befa_spec("raw", "unit_vector", 2L))
  expect_lt(abs(uv$logml), 0.1)
  expect_identical(uv$setup$info$coordinates, "polar")
  # Normal prior (Cartesian coordinates)
  normal <- .prior_logml(.befa_spec("cov", "normal", 2L))
  expect_lt(abs(normal$logml), 0.1)
  expect_identical(normal$setup$info$coordinates, "cartesian")
  # One factor: p(y) = 2 p_+(y)
  uni <- .prior_logml(.befa_spec("cor", "unit_vector", 1L), M = 1L)
  expect_lt(abs(uni$logml), 0.1)
  expect_equal(uni$setup$info$log_vol, log(2))
})

test_that("priors that rotations change are refused", {
  spec <- .befa_spec("cov", "normal", 2L)
  prep <- .prepared_for(spec)
  sd <- .stan_data(spec, prep, 2L, prepare_befa_priors(list(lambda = c(0.5, 1)),
                                                        "cov", "normal", 2L))
  draws <- list(pars = list(Lambda = matrix(1, 4, 10), psi = matrix(1, 4, 5)),
                G = diag(5))
  expect_error(.bridge_setup(spec, sd, draws$pars, 2L, draws$G), "location 0")
  spec <- .befa_spec("cor", "unit_vector", 1L)
  sd <- .stan_data(spec, .prepared_for(spec), 1L,
                   prepare_befa_priors(list(lambda = c(2, 1)), "cor", "unit_vector", 1L))
  expect_error(.bridge_setup(spec, sd, list(lambda_std = matrix(0.5, 4, 5)), 2L, diag(5)),
               "equal shapes")
})

test_that("befa_bridge_sampler() gives the same evidence with any anchors", {
  skip_on_cran()
  skip_if_not_installed("bridgesampling")
  skip_if_not_installed("stanli", "0.25.0")
  skip_if_not(stanli::stanli_available(), "the stanli runtime is not installed")
  withr::local_options(BayesianEFA.cache_dir = withr::local_tempdir())
  fit <- suppressWarnings(befa(.simulated_data(N = 300), 2L, backend = "stanli",
                               chains = 2, iter_warmup = 500, iter_sampling = 1500,
                               seed = 4, verbose = FALSE,
                               compute_fit_indices = FALSE,
                               compute_reliability = FALSE))
  auto <- befa_bridge_sampler(fit)
  other <- befa_bridge_sampler(fit, anchors = c(2L, 4L))
  expect_s3_class(auto, "bridge")
  expect_true(is.finite(auto$logml))
  expect_identical(attr(other, "befa")$anchors, c(2L, 4L))
  expect_lt(abs(auto$logml - other$logml), 0.5)
  expect_error(befa_bridge_sampler(fit, anchors = 1L), "2 distinct")

  old <- fit
  old$fit$draws <- posterior::subset_draws(old$fit$draws, variable = "lp__")
  expect_error(befa_bridge_sampler(old), "does not keep the draws")
  approx <- fit
  approx$method <- "advi"
  expect_warning(befa_bridge_sampler(approx), "not from MCMC")
})

test_that("weak or singular anchors are reported", {
  spec <- .befa_spec("cor", "unit_vector", 2L)
  prep <- .prepared_for(spec)
  sd <- .stan_data(spec, prep, 2L, .prior_for(spec))
  draws <- .prior_draws(spec, sd, S = 40L)
  expect_warning(.bridge_setup(spec, sd, draws$pars, 2L, G = diag(0.01, 5L)),
                 "Weak anchors")
  Z <- matrix(c(1, 1, 0.3, 0.5, 0.5, 0.2), 3, 2)
  expect_error(.rotate_to_T(Z, 1:3, 2L), "singular")
})
