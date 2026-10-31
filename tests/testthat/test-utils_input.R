# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Validate rsp arguments
# ─────────────────────────────────────────────────────────────────────────────

test_that("validate_rsp_args returns defaults when NULL", {
  result <- validate_rsp_args(NULL)
  expect_equal(result$max_iter, 1000)
  expect_equal(result$threshold, 1e-6)
})

test_that("validate_rsp_args merges user values with defaults", {
  result <- validate_rsp_args(list(max_iter = 500))
  expect_equal(result$max_iter, 500)
  expect_equal(result$threshold, 1e-6) # Default preserved
})

test_that("validate_rsp_args errors on non-list input", {
  expect_error(validate_rsp_args("not_a_list"), "must be a named list")
})

test_that("validate_rsp_args warns on unknown arguments", {
  expect_warning(validate_rsp_args(list(unknown_arg = 999)), "Unknown arguments")
})

test_that("validate_rsp_args errors on invalid max_iter", {
  expect_error(validate_rsp_args(list(max_iter = -1)), "positive integer")
  expect_error(validate_rsp_args(list(max_iter = "abc")), "positive integer")
})

test_that("validate_rsp_args errors on invalid threshold", {
  expect_error(validate_rsp_args(list(threshold = 0)), "positive numeric")
  expect_error(validate_rsp_args(list(threshold = -0.001)), "positive numeric")
})

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Validate loo_args
# ─────────────────────────────────────────────────────────────────────────────

test_that("validate_loo_args returns defaults when NULL", {
  result <- validate_loo_args(NULL)
  expect_equal(result$r_eff, FALSE)
  expect_equal(result$cores, 1)
})

test_that("validate_loo_args updates defaults with user input", {
  result <- validate_loo_args(list(r_eff = FALSE, cores = 4))
  expect_equal(result$r_eff, FALSE)
  expect_equal(result$cores, 4)
})

test_that("validate_loo_args warns on unknown arguments", {
  expect_warning(validate_loo_args(list(unknown = "foo")), "Unknown arguments")
})

test_that("validate_loo_args errors on invalid input", {
  expect_error(validate_loo_args(list(r_eff = "yes")), "logical value")
  expect_error(validate_loo_args(list(cores = -1)), "positive integer")
})


# ─────────────────────────────────────────────────────────────────────────────
