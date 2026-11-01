# The C++ log density must give, for every structure, the log density of the
# Stan program (with its constants and Jacobians), which .reference_lp()
# computes row by row in R and test-stan-lp.R checks against Stan.

specs <- .all_specs()
for (label in names(specs)) {
  spec <- specs[[label]]
  if (spec$loadings == "none") next
  test_that(paste("the C++ log density of", label, "matches the reference"), {
    prep <- .prepared_for(spec)
    prior <- .prior_for(spec)
    M <- .spec_factors(spec)
    pars <- .random_pars(spec, prep$J, M)
    item_sd <- if (spec$loadings == "normal") .item_sd(prep)
    u <- drop(.unconstrain(spec, lapply(pars, function(x) matrix(x, nrow = 1))))
    density <- .befa_density(spec, .stan_data(spec, prep, M, prior))
    expect_equal(befa_log_density_cpp(density, u),
                 .reference_lp(spec, prep, pars, prior, item_sd), tolerance = 1e-8)
  })
}

test_that("summary statistics give the same C++ log density as the rows", {
  Y <- .simulated_data()
  spec <- .befa_spec("raw", "unit_vector", 2L)
  prior <- .prior_for(spec)
  pars <- .random_pars(spec, 5L, 2L)
  u <- drop(.unconstrain(spec, lapply(pars, function(x) matrix(x, nrow = 1))))
  rows <- .prepare_befa_data(Y, "raw")
  summ <- .prepare_befa_data(model_type = "raw", sample_nobs = nrow(Y),
                             sample_cov = stats::cov(Y), sample_mean = colMeans(Y))
  lp_rows <- befa_log_density_cpp(.befa_density(spec, .stan_data(spec, rows, 2L, prior)), u)
  lp_summ <- befa_log_density_cpp(.befa_density(spec, .stan_data(spec, summ, 2L, prior)), u)
  expect_equal(lp_summ, lp_rows, tolerance = 1e-8)
})

test_that("the prior-only density drops the likelihood and nothing else", {
  spec <- .befa_spec("cov", "unit_vector", 2L, "multiple")
  prep <- .prepared_for(spec)
  prior <- .prior_for(spec)
  sd <- .stan_data(spec, prep, 2L, prior)
  pars <- .random_pars(spec, prep$J, 2L)
  u <- drop(.unconstrain(spec, lapply(pars, function(x) matrix(x, nrow = 1))))
  full <- befa_log_density_cpp(.befa_density(spec, sd), u)
  prior_only <- befa_log_density_cpp(.befa_density(spec, sd, prior_only = TRUE), u)
  # The likelihood, row by row, at the same parameters
  L <- pars$Z * sqrt(pars$h2 / rowSums(pars$Z^2))
  Sigma <- (tcrossprod(L) + diag(1 - pars$h2)) * tcrossprod(pars$sigma)
  ll <- sum(vapply(seq_len(nrow(prep$Y)), function(i) {
    o <- which(!is.na(prep$Y[i, ]))
    mvnfast::dmvn(prep$Y[i, o], rep(0, length(o)), Sigma[o, o, drop = FALSE], log = TRUE)
  }, 0))
  expect_equal(full - prior_only, ll, tolerance = 1e-8)
})
