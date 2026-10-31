test_that("the ML EFA reproduces stats::factanal()", {
  Y <- .simulated_data(N = 300)
  R <- stats::cor(Y)
  efa <- .ml_efa(R, 2L)
  fa <- stats::factanal(covmat = R, factors = 2, n.obs = 300, rotation = "none")
  expect_true(efa$converged)
  expect_equal(efa$uniquenesses, unname(fa$uniquenesses), tolerance = 1e-4)
  # Loadings agree up to rotation, so compare their cross-products.
  expect_equal(tcrossprod(efa$loadings), unname(tcrossprod(unclass(fa$loadings))),
               tolerance = 1e-3)
})

test_that("a failed ML EFA falls back to principal axes with SMCs", {
  R <- stats::cor(.simulated_data())
  local_mocked_bindings(.ml_efa = function(R, M) stop("no convergence"))
  sol <- .start_solution(R, 2L)
  expect_false(sol$converged)
  expect_identical(dim(sol$loadings), c(5L, 2L))
  expect_equal(sol$uniquenesses, 1 - rowSums(sol$loadings^2))
})

test_that("the starting moments come from rows, the EM or summary statistics", {
  Y <- .simulated_data()
  m <- .start_moments(.prepare_befa_data(Y, "raw"))
  expect_equal(m$mean, colMeans(Y))
  expect_equal(m$sd, apply(Y, 2, .ml_sd))
  prep <- .prepare_befa_data(.simulated_data(missing = TRUE), "raw", "FIML")
  em <- estimate_h1_moments_em(prep$Y)
  expect_equal(.start_moments(prep)$mean, em$Mu)
  s <- .prepare_befa_data(model_type = "raw", sample_nobs = 40, sample_cov = stats::cov(Y),
                          sample_mean = colMeans(Y))
  expect_equal(.start_moments(s)$mean, colMeans(Y))
  sc <- .prepare_befa_data(model_type = "cor", sample_nobs = 40, sample_cor = stats::cor(Y))
  expect_equal(.start_moments(sc)$mean, rep(0, 5))
})

test_that("starting values fit the parameters of every structure", {
  for (spec in .spec_grid()) {
    prep <- .prepared_for(spec)
    prior <- .prior_for(spec)
    M <- .spec_factors(spec)
    inits <- .start_values(spec, prep, M, prior, chains = 2L, seed = 9L)
    expect_length(inits, 2L)
    expect_setequal(names(inits[[1L]]), names(.random_pars(spec, prep$J, M)))
    expect_false(identical(inits[[1L]], inits[[2L]]))
  }
})

test_that("unit-vector starts put the norms of Z at the mode of the repulsion", {
  spec <- .befa_spec("cor", "unit_vector", 2L)
  prior <- .prior_for(spec)
  init <- .start_values(spec, .prepared_for(spec), 2L, prior, 1L, 3L)[[1L]]
  expect_equal(sqrt(rowSums(init$Z^2)), rep(sqrt(prior$xi + 1), 5))
  expect_true(all(init$h2 >= 0.05 & init$h2 <= 0.95))
})

test_that("starting values are reproducible and leave the random stream alone", {
  spec <- .befa_spec("raw", "normal", 2L)
  prep <- .prepared_for(spec)
  prior <- .prior_for(spec)
  set.seed(1)
  before <- .Random.seed
  a <- .start_values(spec, prep, 2L, prior, 2L, 5L)
  expect_identical(.Random.seed, before)
  expect_identical(a, .start_values(spec, prep, 2L, prior, 2L, 5L))
  expect_true(all(a[[1L]]$psi > 0))
})
