test_that("the C++ EM gives the v0.1.0 saturated moments", {
  fx <- .read_fixture("em_h1_moments")
  res <- estimate_h1_moments_em(fx$inputs$Y)
  # Grouping rows by pattern only changes the order of the sums.
  expect_equal(res$Mu, fx$output$Mu, tolerance = 1e-8)
  expect_equal(res$Sigma, fx$output$Sigma, tolerance = 1e-8)
})

test_that("complete data give the ML moments after one EM step", {
  set.seed(31)
  Y <- matrix(stats::rnorm(60 * 3), ncol = 3)
  res <- em_h1_moments_cpp(Y)
  expect_equal(res$Mu, colMeans(Y), tolerance = 1e-12)
  expect_equal(res$Sigma, stats::cov(Y) * 59 / 60, tolerance = 1e-12)
  # The second step changes nothing, so the loop stops there.
  expect_identical(res$iterations, 2L)
})

test_that("the EM stops at max_iter", {
  set.seed(32)
  Y <- matrix(stats::rnorm(40 * 3), ncol = 3)
  Y[sample(length(Y), 20)] <- NA
  expect_identical(em_h1_moments_cpp(Y, max_iter = 1L)$iterations, 1L)
})

test_that("rows and items with no observed value keep the moments finite", {
  set.seed(33)
  Y <- matrix(stats::rnorm(30 * 3), ncol = 3)
  Y[5, ] <- NA
  Y[, 3] <- NA
  res <- estimate_h1_moments_em(Y)
  expect_true(all(is.finite(res$Mu)))
  expect_true(all(is.finite(res$Sigma)))
})

test_that("a singular covariance gets the small ridge", {
  set.seed(34)
  x <- stats::rnorm(50)
  Y <- cbind(x, x, stats::rnorm(50))
  res <- estimate_h1_moments_em(Y)
  expect_gt(min(eigen(res$Sigma, symmetric = TRUE)$values), 0)
})
