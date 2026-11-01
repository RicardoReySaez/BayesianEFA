test_that("the extra variables of every backend get the same names", {
  vars <- c("a", "lp__", "log_p__", "log_g__", "path__", "treedepth__")
  draws <- posterior::as_draws_array(array(seq_len(12), c(2, 1, 6),
                                           dimnames = list(NULL, NULL, vars)))
  out <- .standardize_draws(draws)
  expect_identical(posterior::variables(out), c("a", "lp__", "lp_approx__"))
  # rstan's lp__ of ADVI is replaced by log_p__.
  expect_equal(as.numeric(out[, , "lp__"]), 5:6)
  plain <- posterior::as_draws_array(array(1:4, c(2, 1, 2),
                                           dimnames = list(NULL, NULL, c("b", "lp__"))))
  expect_identical(posterior::variables(.standardize_draws(plain)), c("b", "lp__"))
})

test_that("durations, iterations and the bar style are read and written", {
  expect_identical(.format_elapsed(12.34), "12.3s")
  expect_identical(.format_elapsed(125), "2m 5s")
  out <- withr::local_tempfile()
  writeLines(c("Chain 1 Iteration:   1 / 40", "noise", "Chain 1 Iteration:  20 / 40"), out)
  job <- list(files = list(out = out, err = tempfile()))
  expect_identical(.reported_iteration(job), 20L)
  expect_identical(.read_lines(tempfile()), character(0))
  writeLines("nothing", out)
  expect_identical(.reported_iteration(job), 0L)
  withr::local_options(cli.unicode = TRUE)
  expect_named(.bar_style(), c("complete", "current", "incomplete"))
  withr::local_options(cli.unicode = FALSE)
  expect_null(.bar_style())
})

test_that("a background job returns its value or its error", {
  skip_on_cran()
  ok <- .launch_job(function(x) x + 1, list(x = 1))
  ok$process$wait()
  expect_identical(.job_value(ok, "Job"), 2)
  .stop_job(ok)
  expect_false(any(file.exists(unlist(ok$files))))

  bad <- .launch_job(function() stop("boom"), list())
  bad$process$wait()
  expect_error(.job_value(bad, "Chain 3"), "Chain 3 failed: boom")

  dead <- .launch_job(function() {
    cat("last words\n")
    quit(save = "no", status = 1)
  }, list())
  dead$process$wait()
  expect_error(.job_value(dead, "Chain 1"), "stopped without a result.*last words")
})

# This function returns a fake MCMC job that reports its iterations as Stan
# does and returns its number.
.fake_chain <- function(i, iters = 4L) {
  list(func = function(i, iters) {
    for (k in seq_len(iters)) {
      cat(sprintf("Chain %d Iteration: %d / %d\n", i, k, iters))
      Sys.sleep(0.05)
    }
    i
  }, args = list(i = i, iters = iters))
}

test_that("jobs run a few at a time and report once they finish", {
  skip_on_cran()
  jobs <- lapply(1:3, .fake_chain)
  expect_message(res <- .run_jobs(jobs, cores = 2L, per_job = 4L, warmup = 2L,
                                  label = "Chain", done = "Sampled 3 chains"),
                 "Sampled 3 chains in")
  expect_identical(res, as.list(1:3))
  expect_silent(.run_jobs(jobs[1], per_job = 4L, quiet = TRUE))
})

test_that("a console that redraws in place shows the bar and the spinner", {
  skip_on_cran()
  withr::local_options(cli.dynamic = TRUE, cli.unicode = TRUE)
  expect_no_error(suppressMessages(
    .run_jobs(lapply(1:2, .fake_chain), per_job = 4L, warmup = 2L, done = "Sampled")
  ))
  expect_no_error(suppressMessages(
    .run_jobs(list(.fake_chain(1L)), label = "ADVI", done = "Ran ADVI")
  ))
})

test_that("a failing job stops the others", {
  skip_on_cran()
  jobs <- list(.fake_chain(1L, iters = 40L),
               list(func = function() stop("bad data"), args = list()))
  expect_error(.run_jobs(jobs, per_job = 40L, label = "Chain", quiet = TRUE),
               "Chain 2 failed: bad data")
})

test_that("the header describes how each method runs", {
  header <- function(...) paste(capture_messages(.stan_header(...)), collapse = "")
  mcmc <- .sampler_settings("mcmc", list(chains = 2, parallel_chains = 4, seed = 1))
  expect_match(header("cmdstanr", mcmc), "cmdstanr")
  expect_match(header("cmdstanr", mcmc),
               "2 chains of 1000 draws after 1000 warmup iterations, 2 at a time")
  advi <- .sampler_settings("advi", list(seed = 1))
  expect_match(header("stanr", advi), "ADVI \\(meanfield\\) with up to 10000")
  path <- .sampler_settings("pathfinder", list(seed = 1))
  expect_match(header("stanli", path), "Pathfinder with 1 path and 1000 draws")
  expect_match(header("cmdstanr", path), "Pathfinder with 4 paths")
})

test_that("chains are joined into one set of draws", {
  one <- function(v) {
    list(draws = posterior::as_draws_array(array(v, c(2, 1, 2),
                                                 dimnames = list(NULL, NULL, c("a", "lp__")))),
         sampler = posterior::as_draws_array(array(0, c(2, 1, 1),
                                                   dimnames = list(NULL, NULL, "divergent__"))),
         seconds = 1)
  }
  s <- .sampler_settings("mcmc", list(seed = 1))
  fit <- .combine_jobs(list(one(1:4), one(5:8)), "stanli", s)
  expect_s3_class(fit, "befa_draws")
  expect_identical(posterior::nchains(fit$draws), 2L)
  expect_identical(posterior::nchains(fit$sampler), 2L)
  expect_identical(fit$seconds, c(1, 1))
  single <- .combine_jobs(list(list(draws = one(1:4)$draws, sampler = NULL, seconds = 2)),
                          "stanli", .sampler_settings("advi", list(seed = 1)))
  expect_null(single$sampler)
  expect_identical(single$method, "advi")
})

# A small two-factor problem shared by the tests that run Stan.
.small_problem <- function() {
  spec <- .befa_spec("cor", "unit_vector", 2L)
  prep <- .prepare_befa_data(.simulated_data(N = 150), "cor")
  prior <- .prior_for(spec)
  list(spec = spec, prep = prep, prior = prior, code = .stan_code(spec),
       data = .stan_data(spec, prep, 2L, prior))
}

test_that("a stanli job samples, approximates and runs Pathfinder in this session", {
  skip_on_cran()
  skip_if_not_installed("stanli", "0.25.0")
  skip_if_not(stanli::stanli_available(), "the stanli runtime is not installed")
  p <- .small_problem()
  file <- withr::local_tempfile(fileext = ".stan")
  writeLines(p$code, file)
  for (method in .befa_methods) {
    s <- .sampler_settings(method, list(seed = 3, iter_warmup = 100, iter_sampling = 50,
                                        refresh = 0)[if (method == "mcmc") 1:4 else c(1, 4)])
    n <- switch(method, mcmc = 1L, advi = 1L, pathfinder = 2L)
    init <- .start_values(p$spec, p$prep, 2L, p$prior, n, 3L)
    if (method != "pathfinder") init <- init[[1L]]
    out <- .stan_job("stanli", method, file, p$data, s, init, 1L)
    expect_true("lp__" %in% posterior::variables(out$draws), label = method)
    expect_identical(!is.null(out$sampler), method == "mcmc")
  }
})

test_that("stanli runs every method through the background processes", {
  skip_on_cran()
  skip_if_not_installed("stanli", "0.25.0")
  skip_if_not(stanli::stanli_available(), "the stanli runtime is not installed")
  withr::local_options(BayesianEFA.cache_dir = withr::local_tempdir())
  p <- .small_problem()
  s <- .sampler_settings("mcmc", list(chains = 2, iter_warmup = 150, iter_sampling = 100,
                                      seed = 11))
  inits <- .start_values(p$spec, p$prep, 2L, p$prior, 2L, s$seed)
  msgs <- capture_messages(fit <- .run_stan(p$code, p$data, inits, "stanli", s))
  expect_match(paste(msgs, collapse = ""), "Sampled 2 chains in")
  expect_identical(unname(dim(fit$draws)[1:2]), c(100L, 2L))
  expect_true(all(posterior::variables(fit$sampler) %in%
                    c("accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__",
                      "divergent__", "energy__")))

  advi <- .sampler_settings("advi", list(seed = 11, draws = 200))
  fit <- .run_stan(p$code, p$data, inits[[1L]], "stanli", advi, quiet = TRUE)
  expect_true(all(c("lp__", "lp_approx__") %in% posterior::variables(fit$draws)))
  expect_identical(unname(posterior::ndraws(fit$draws)), 200L)
})

test_that("cmdstanr and stanr fit the same model (opt-in, they compile)", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("BAYESIANEFA_TEST_BACKENDS"), "true"),
              "set BAYESIANEFA_TEST_BACKENDS=true to compile with cmdstanr and stanr")
  withr::local_options(BayesianEFA.cache_dir = withr::local_tempdir())
  p <- .small_problem()
  for (backend in c("cmdstanr", "stanr")) {
    if (!.backend_status()[[backend]]$available) next
    for (method in .befa_methods) {
      s <- .sampler_settings(method, list(seed = 5))
      n <- switch(method, mcmc = s$chains, advi = 1L, pathfinder = s$num_paths)
      inits <- .start_values(p$spec, p$prep, 2L, p$prior, n, s$seed)
      if (method == "advi") inits <- inits[[1L]]
      fit <- .run_stan(p$code, p$data, inits, backend, s, quiet = TRUE)
      h2 <- colMeans(posterior::as_draws_matrix(posterior::subset_draws(fit$draws, "h2")))
      expect_true(all(h2 > 0.05 & h2 < 0.95), label = paste(backend, method))
    }
  }
})

# This function returns a Pathfinder-like fit: one chain of candidate draws
# from N(0, 1) with their target density and that of a wider normal.
.fake_pathfinder <- function(n = 4000, seed = 31) {
  set.seed(seed)
  x <- stats::rnorm(n, 0, 1.2)
  draws <- posterior::as_draws_array(array(
    c(x, stats::dnorm(x, log = TRUE), stats::dnorm(x, 0, 1.2, log = TRUE)),
    dim = c(n, 1, 3), dimnames = list(NULL, NULL, c("mu", "lp__", "lp_approx__"))))
  structure(list(draws = draws, method = "pathfinder"), class = "befa_draws")
}

test_that("multi-path Pathfinder is resampled after its Pareto k", {
  settings <- list(psis_resample = TRUE, draws = 1000L, seed = 2L)
  fit <- .pathfinder_resample(.fake_pathfinder(), settings)
  expect_lt(fit$pareto_k, 0.5)
  expect_identical(posterior::ndraws(fit$draws), 1000L)
  # The PSIS weights move the draws from N(0, 1.2^2) towards N(0, 1).
  expect_equal(stats::sd(.stacked_draws(fit$draws, "mu")), 1, tolerance = 0.08)
  expect_identical(.pathfinder_resample(.fake_pathfinder(), settings)$draws,
                   fit$draws)

  settings$psis_resample <- FALSE
  kept <- .pathfinder_resample(.fake_pathfinder(), settings)
  expect_identical(posterior::ndraws(kept$draws), 4000L)
  expect_identical(kept$pareto_k, fit$pareto_k)
})

test_that("Pathfinder draws without both densities are left as they are", {
  fit <- .fake_pathfinder()
  fit$draws <- posterior::subset_draws(fit$draws, variable = c("mu", "lp__"))
  settings <- list(psis_resample = TRUE, draws = 1000L, seed = 2L)
  expect_identical(.pathfinder_resample(fit, settings), fit)
  fit <- .fake_pathfinder()
  fit$draws[1, 1, "lp__"] <- -Inf
  expect_null(.pathfinder_resample(fit, settings)$pareto_k)
})

test_that(".run_stan() resamples the Pathfinder of the compiled backends", {
  candidates <- .fake_pathfinder()
  local_mocked_bindings(
    .compiled_model = function(...) list(stan_file = "model.stan"),
    .run_jobs = function(...) list(list(draws = candidates$draws, sampler = NULL,
                                        seconds = 1))
  )
  s <- .sampler_settings("pathfinder", list(seed = 3, refresh = 0))
  fit <- .run_stan("code", list(), list(list()), "cmdstanr", s, quiet = TRUE)
  expect_identical(posterior::ndraws(fit$draws), 1000L)
  expect_lt(fit$pareto_k, 0.5)
  kept <- .run_stan("code", list(), list(list()), "stanli", s, quiet = TRUE)
  expect_identical(posterior::ndraws(kept$draws), 4000L)
  expect_null(kept$pareto_k)
})
