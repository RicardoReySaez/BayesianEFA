# This function returns a fake MCMC fit whose draws are well mixed, with
# the given number of divergent transitions and maximum-depth hits.
.fake_mcmc <- function(divergent = 0L, deep = 0L, n_iter = 500L, seed = 71) {
  set.seed(seed)
  draws <- posterior::as_draws_array(array(stats::rnorm(n_iter * 4 * 2), c(n_iter, 4, 2),
                                           dimnames = list(NULL, NULL, c("h2[1]", "h2[2]"))))
  div <- rep(0, n_iter * 4)
  div[seq_len(divergent)] <- 1
  depth <- rep(5, n_iter * 4)
  depth[seq_len(deep)] <- 15
  sampler <- posterior::as_draws_array(array(c(div, depth), c(n_iter, 4, 2),
                                             dimnames = list(NULL, NULL, c("divergent__", "treedepth__"))))
  list(draws = draws, fit = list(method = "mcmc", backend = "stanli", sampler = sampler))
}

test_that("well mixed chains pass every check", {
  f <- .fake_mcmc()
  d <- .befa_diagnostics(f$draws, f$fit, list(max_treedepth = 15L))
  expect_identical(c(d$divergent, d$treedepth), c(0L, 0L))
  expect_lt(d$max_rhat, 1.01)
  expect_gt(min(d$min_ess_bulk, d$min_ess_tail), 400)
  expect_length(d$issues, 0L)
  expect_false(.befa_warn_issues(d))
  msgs <- paste(capture_messages(.befa_report(d, list(max_treedepth = 15L))), collapse = "")
  expect_match(msgs, "No divergent transitions")
  expect_match(msgs, "Largest R-hat .* below 1.01")
  expect_match(msgs, "Smallest effective sample size .* above 400")
})

test_that("sampler problems are named, reported and turned into a warning", {
  f <- .fake_mcmc(divergent = 3L, deep = 2L, n_iter = 50L)
  d <- .befa_diagnostics(f$draws, f$fit, list(max_treedepth = 15L))
  expect_identical(c(d$divergent, d$treedepth), c(3L, 2L))
  expect_match(d$issues[1], "3 divergent transitions")
  expect_match(d$issues[2], "2 transitions hit the maximum tree depth")
  expect_true(any(grepl("effective sample size", d$issues)))
  msgs <- paste(capture_messages(.befa_report(d, list(max_treedepth = 15L))), collapse = "")
  expect_match(msgs, "3 divergent transitions")
  expect_match(msgs, "2 transitions reached the maximum tree depth \\(15\\)")
  expect_match(msgs, "below 400")
  expect_warning(.befa_warn_issues(d), "adapt_delta.*max_treedepth.*iter_sampling")
})

test_that("chains that disagree raise R-hat", {
  f <- .fake_mcmc()
  arr <- unclass(f$draws)
  arr[, 1, ] <- arr[, 1, ] + 5
  d <- .befa_diagnostics(posterior::as_draws_array(arr), f$fit, list(max_treedepth = 15L))
  expect_gt(d$max_rhat, 1.01)
  expect_match(paste(capture_messages(.befa_report(d, list(max_treedepth = 15L))),
                     collapse = ""), "above 1.01")
})

test_that("a fit without sampler diagnostics still reports R-hat and ESS", {
  f <- .fake_mcmc()
  f$fit$sampler <- NULL
  d <- .befa_diagnostics(f$draws, f$fit, list(max_treedepth = 15L))
  expect_true(is.na(d$divergent))
  msgs <- paste(capture_messages(.befa_report(d, list(max_treedepth = 15L))), collapse = "")
  expect_false(grepl("divergent", msgs))
})

# This function returns fake draws of an approximation whose log densities
# differ from the target's by noise of the given size.
.fake_approx <- function(noise, n = 1000L) {
  set.seed(72)
  lp <- stats::rnorm(n)
  arr <- array(c(lp, lp + stats::rnorm(n, 0, noise)), c(n, 1, 2),
               dimnames = list(NULL, NULL, c("lp__", "lp_approx__")))
  posterior::as_draws_array(arr)
}

test_that("approximations are judged by the Pareto k of their importance ratios", {
  good <- .befa_diagnostics(NULL, list(method = "advi", backend = "cmdstanr",
                                       draws = .fake_approx(0.1)), list())
  expect_lt(good$pareto_k, 0.7)
  expect_length(good$issues, 0L)
  expect_match(paste(capture_messages(.befa_report(good, list())), collapse = ""), "below 0.7")

  bad <- .befa_diagnostics(NULL, list(method = "pathfinder", backend = "cmdstanr",
                                      draws = .fake_approx(5)), list())
  expect_gt(bad$pareto_k, 0.7)
  expect_match(paste(capture_messages(.befa_report(bad, list())), collapse = ""), "above 0.7")
  expect_warning(.befa_warn_issues(bad), "method = \"mcmc\"")

  # The Pareto k of all the candidates of Pathfinder, before resampling, wins.
  resampled <- .befa_diagnostics(NULL, list(method = "pathfinder", backend = "cmdstanr",
                                            draws = .fake_approx(5), pareto_k = 0.3),
                                 list())
  expect_identical(resampled$pareto_k, 0.3)
})

test_that("the Pareto k is missing when the backend lacks a density", {
  only_approx <- posterior::as_draws_array(array(0, c(10, 1, 1),
                                                 dimnames = list(NULL, NULL, "lp_approx__")))
  d <- .befa_diagnostics(NULL, list(method = "advi", backend = "stanr",
                                    draws = only_approx), list())
  expect_true(is.na(d$pareto_k))
  expect_match(paste(capture_messages(.befa_report(d, list())), collapse = ""),
               "does not report the density")
  infinite <- .fake_approx(0.1)
  infinite[1, 1, "lp__"] <- -Inf
  expect_true(is.na(.approximation_k(infinite)))
})
