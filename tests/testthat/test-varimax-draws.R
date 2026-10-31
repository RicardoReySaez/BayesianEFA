test_that("the C++ varimax reproduces stats::varimax draw by draw", {
  set.seed(11)
  for (M in 2:5) {
    J <- 3L * M
    draws <- matrix(stats::rnorm(40 * J * M, 0, 0.5), nrow = 40)
    rotated <- varimax_draws_cpp(draws, J, M)
    expected <- t(apply(draws, 1, function(x) {
      c(stats::varimax(matrix(x, ncol = M), normalize = FALSE)$loadings)
    }))
    # Both call LAPACK's dgesdd and stop at the same iteration, so only
    # rounding separates them.
    expect_equal(rotated, expected, tolerance = 1e-10)
  }
})

test_that("a different eps changes the stopping point as in stats::varimax", {
  set.seed(12)
  draws <- matrix(stats::rnorm(20 * 12), nrow = 20)
  rotated <- varimax_draws_cpp(draws, 6L, 2L, eps = 1e-2)
  expected <- t(apply(draws, 1, function(x) {
    c(stats::varimax(matrix(x, ncol = 2), normalize = FALSE, eps = 1e-2)$loadings)
  }))
  expect_equal(rotated, expected, tolerance = 1e-10)
})

test_that("a single factor is returned unchanged", {
  draws <- matrix(seq_len(10) / 10, nrow = 2)
  expect_identical(varimax_draws_cpp(draws, 5L, 1L), draws)
})

test_that("draws of the wrong width are refused", {
  draws <- matrix(0, nrow = 3, ncol = 7)
  expect_error(varimax_draws_cpp(draws, 3L, 2L), "J \\* M = 6")
})
