test_that("each method gets its defaults", {
  mcmc <- .sampler_settings("mcmc", list(seed = 1))
  expect_identical(mcmc[c("chains", "iter_warmup", "iter_sampling", "max_treedepth")],
                   list(chains = 4L, iter_warmup = 1000L, iter_sampling = 1000L,
                        max_treedepth = 15L))
  expect_identical(mcmc$adapt_delta, 0.9)
  expect_identical(mcmc$refresh, 20L)
  expect_identical(mcmc$method, "mcmc")

  advi <- .sampler_settings("advi", list(seed = 1))
  expect_identical(advi$algorithm, "meanfield")
  expect_identical(advi$refresh, 100L)
  path <- .sampler_settings("pathfinder", list(seed = 1))
  expect_identical(path$num_paths, 4L)
  expect_true(path$psis_resample)
})

test_that("the rstan names of the MCMC options are understood", {
  s <- .sampler_settings("mcmc", list(iter = 3000, warmup = 1000, cores = 2,
                                      control = list(adapt_delta = 0.9), seed = 1))
  expect_identical(c(s$iter_warmup, s$iter_sampling, s$parallel_chains), c(1000L, 2000L, 2L))
  expect_identical(s$adapt_delta, 0.9)
  # Without warmup, rstan keeps half of iter for it.
  half <- .sampler_settings("mcmc", list(iter = 400, seed = 1))
  expect_identical(c(half$iter_warmup, half$iter_sampling), c(200L, 200L))
  expect_identical(.sampler_settings("mcmc", list(warmup = 50, seed = 1))$iter_warmup, 50L)
  expect_error(.sampler_settings("mcmc", list(iter = 100, warmup = 100)), "must exceed")
  expect_error(.sampler_settings("mcmc", list(control = list(stepsize = 1))),
               "accepts `adapt_delta` and `max_treedepth`, not `stepsize`")
})

test_that("unknown or invalid options are errors", {
  expect_error(.sampler_settings("mcmc", list(chainz = 4)), "`chainz` is not an option of MCMC")
  expect_error(.sampler_settings("mcmc", list(chains = 0)), "positive whole number")
  expect_error(.sampler_settings("mcmc", list(adapt_delta = 1)), "between 0 and 1")
  expect_error(.sampler_settings("advi", list(algorithm = "laplace")), "must be one of")
  expect_error(.sampler_settings("pathfinder", list(num_paths = -1)), "positive whole number")
  expect_error(.sampler_settings("advi", list(refresh = -1)), "`refresh` must be zero")
  expect_error(.sampler_settings("mcmc", list(seed = -2)), "non-negative whole number")
})

test_that("a seed is drawn when none is given, and zero is a valid seed", {
  set.seed(1)
  s <- .sampler_settings("mcmc")
  expect_type(s$seed, "integer")
  expect_identical(.sampler_settings("advi", list(seed = 0))$seed, 0L)
})

test_that("R CMD check limits the default number of parallel chains", {
  withr::local_envvar(`_R_CHECK_LIMIT_CORES_` = "TRUE")
  expect_lte(.default_cores(), 2L)
  withr::local_envvar(`_R_CHECK_LIMIT_CORES_` = "false")
  expect_gte(.default_cores(), 1L)
})
