# Omega coefficients from the standardized draws

test_that("omega follows the standardized loadings draw by draw", {
  object <- .fake_befa()
  rel <- befa_reliability(object)
  expect_s3_class(rel, "befa_reliability")
  expect_identical(posterior::variables(rel$omega_draws),
                   c("omega_total", "omega_F1", "omega_F2"))
  expect_identical(unname(dim(rel$omega_draws)), c(100L, 2L, 3L))

  std <- .metric_draws(object, "std")
  L <- matrix(.stacked_draws(std, "Lambda_std")[7, ], 5, 2)
  psi <- sum(1 - rowSums(L^2))
  sums <- colSums(L)
  omega <- .stacked_draws(rel$omega_draws, posterior::variables(rel$omega_draws))[7, ]
  expect_equal(omega, c(sum(sums^2) / (sum(sums^2) + psi), sums^2 / (sums^2 + psi)))
  expect_true(all(omega > 0 & omega < 1))
})

test_that("omega is the same in every metric of the same solution", {
  cor <- befa_reliability(.fake_befa("cor"))
  cov <- befa_reliability(.fake_befa("cov", lambda_prior = "normal"))
  raw <- befa_reliability(.fake_befa("raw"))
  expect_equal(unclass(cov$omega_draws), unclass(cor$omega_draws))
  expect_equal(unclass(raw$omega_draws), unclass(cor$omega_draws))
})

test_that("the summaries hold the mean, SD and credible interval", {
  rel <- befa_reliability(.fake_befa(), probs = c(0.1, 0.9))
  total <- as.vector(.stacked_draws(rel$omega_draws, "omega_total"))
  expect_named(rel$total, c("Estimate", "SD", "Lower", "Upper"))
  expect_equal(rel$total$Estimate, mean(total))
  expect_equal(rel$total$SD, stats::sd(total))
  expect_equal(c(rel$total$Lower, rel$total$Upper),
               stats::quantile(total, c(0.1, 0.9), names = FALSE))
  expect_identical(rownames(rel$subscales), c("F1", "F2"))
  expect_true(all(rel$subscales$Lower < rel$subscales$Estimate &
                    rel$subscales$Estimate < rel$subscales$Upper))
  expect_identical(rel$probs, c(0.1, 0.9))
})

test_that("stored coefficients are returned and probabilities are checked", {
  object <- .fake_befa()
  object$reliability <- "stored"
  expect_identical(befa_reliability(object), "stored")
  object$reliability <- NULL
  expect_error(befa_reliability(object, probs = 0.5), "two probabilities")
  expect_error(befa_reliability(object, probs = c(0, 1)), "two probabilities")
})
