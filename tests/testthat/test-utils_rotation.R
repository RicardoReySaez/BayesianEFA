# ═══════════════════════════════════════════════════════════════════════════════
# ═══════════════════════════════════════════════════════════════════════════════
#                    SECTION 2: rsp_align with varimax
# ═══════════════════════════════════════════════════════════════════════════════

# --- Setup: create synthetic factor draws ---
setup_factor_draws <- function(n_draws = 100, n_items = 6, n_factors = 2) {
  set.seed(42)

  # True loading matrix
  true_lambda <- matrix(c(
    0.8, 0.1,
    0.75, 0.1,
    0.7, 0.1,
    0.1, 0.8,
    0.1, 0.75,
    0.1, 0.7
  ), nrow = 6, ncol = 2, byrow = TRUE)

  # Generate draws with small noise
  draws <- matrix(NA, nrow = n_draws, ncol = n_items * n_factors)
  for (s in 1:n_draws) {
    noise <- matrix(rnorm(n_items * n_factors, 0, 0.1), nrow = n_items, ncol = n_factors)
    draws[s, ] <- c(true_lambda + noise)
  }

  list(
    draws = draws,
    true_lambda = true_lambda,
    n_items = n_items,
    n_factors = n_factors
  )
}

test_that("rsp_align: works with varimax rotation", {
  test_data <- setup_factor_draws()

  result <- rsp_align(
    lambda_draws = test_data$draws,
    n_items = test_data$n_items,
    n_chains = 1,
    format = "column_major",
    n_factors = test_data$n_factors
  )

  expect_type(result, "list")
  expect_true("Lambda_hat_mcmc" %in% names(result))
  expect_true("objective" %in% names(result))
  expect_equal(dim(result$Lambda_hat_mcmc), dim(test_data$draws))
})

test_that("rsp_align: M=1 skips rotation", {
  set.seed(42)

  n_draws <- 50
  n_items <- 5
  draws <- matrix(rnorm(n_draws * n_items), nrow = n_draws)

  result <- suppressWarnings(rsp_align(
    lambda_draws = draws,
    n_chains = 1,
    format = "column_major",
    n_items = n_items,
    n_factors = 1
  ))

  expect_type(result, "list")
  expect_equal(ncol(result$Lambda_hat_mcmc), n_items)
})

# ─────────────────────────────────────────────────────────────────────────────
