# Factor scores drawn from the final draws of a fit

test_that("the scores are draws of every person on every factor", {
  object <- .fake_befa()
  set.seed(1)
  scores <- befa_factor_scores(object)
  expect_s3_class(scores, "befa_factor_scores")
  expect_s3_class(scores$eta_draws, "draws_array")
  expect_identical(unname(dim(scores$eta_draws)), c(100L, 2L, 120L))
  expect_identical(posterior::variables(scores$eta_draws)[c(1, 60, 61)],
                   c("eta[1,1]", "eta[60,1]", "eta[1,2]"))
  expect_identical(scores$summary$variable, posterior::variables(scores$eta_draws))
  expect_null(befa_factor_scores(object, post_summaries = FALSE)$summary)
})

test_that("the scores follow the seed and match their conditional mean", {
  object <- .fake_befa()
  set.seed(3)
  first <- befa_factor_scores(object, FALSE)$eta_draws
  set.seed(3)
  expect_identical(befa_factor_scores(object, FALSE)$eta_draws, first)

  # Draw by draw, the scores scatter around Sigma_eta L' Psi^-1 y; their
  # average over many draws of the same parameters is close to it.
  p <- .fitted_parameters(object)
  L <- matrix(p$lambda[1, ], 5, 2)
  psi <- p$psi[1, ]
  Sigma_eta <- solve(diag(2) + crossprod(L / psi, L))
  mu <- object$data$Y %*% (L / psi) %*% Sigma_eta
  S <- 4000
  set.seed(4)
  eta <- compute_factor_scores_cpp(
    matrix(p$lambda[1, ], S, 10, byrow = TRUE), matrix(psi, S, 5, byrow = TRUE),
    matrix(0, S, 5), object$data$Y, 5L, 2L, FALSE, list())
  expect_lt(max(abs(matrix(colMeans(eta), ncol = 2) - mu)), 0.08)
  expect_equal(stats::cov(eta[, c(1, 61)]), Sigma_eta, tolerance = 0.1)
})

test_that("the scores of a person condition on the items they observe", {
  object <- .fake_befa("raw", missing = "FIML")
  set.seed(5)
  eta <- .factor_score_draws(object)
  expect_true(all(is.finite(eta)))
  expect_identical(posterior::nvariables(eta), 2L * nrow(object$data$Y))
})

test_that("stored scores are returned and summary statistics have none", {
  object <- .fake_befa()
  object$eta <- befa_factor_scores(object, FALSE)$eta_draws
  expect_identical(befa_factor_scores(object)$eta_draws, object$eta)
  expect_true("eta[1,1]" %in% posterior::variables(.metric_draws(object, "std")))

  summ <- .fake_befa()
  summ$data$summary <- TRUE
  expect_error(befa_factor_scores(summ), "need the raw `data`")
  expect_error(befa_factor_scores(list()), "fit returned by befa")
})

test_that("a person without observed items gets draws from the prior", {
  # RcppArmadillo draws from its own generator, seeded by set.seed(), so the
  # draws are reproducible but are not those of stats::rnorm().
  S <- 5000
  draw <- function() {
    set.seed(9)
    compute_factor_scores_cpp(matrix(0.5, S, 4), matrix(0.5, S, 2),
                              matrix(0, S, 2), matrix(NA_real_, 1, 2), 2L, 2L,
                              TRUE, list(integer(0)))
  }
  eta <- draw()
  expect_identical(draw(), eta)
  expect_lt(max(abs(colMeans(eta))), 0.05)
  expect_lt(max(abs(stats::cov(eta) - diag(2))), 0.06)
})
