specs <- .all_specs()

for (label in names(specs)) {
  test_that(paste("the Stan code of", label, "is stable"), {
    expect_snapshot(cat(.stan_code(specs[[label]])))
  })
}

test_that("the code depends on the structure alone", {
  a <- .stan_code(.befa_spec("cor", "unit_vector", 2L))
  b <- .stan_code(.befa_spec("cor", "unit_vector", 5L))
  expect_identical(a, b)
  expect_false(identical(a, .stan_code(.befa_spec("cov", "unit_vector", 2L))))
})

test_that("only the structures that need them get the optional blocks", {
  uni <- .stan_code(.befa_spec("cor", "unit_vector", 1L))
  expect_match(uni, "functions {", fixed = TRUE)
  uv <- .stan_code(specs[["cor-uv-single-cholesky"]])
  expect_false(grepl("functions {", uv, fixed = TRUE))
  expect_match(.stan_code(.befa_null_spec("cov", "multiple")), "root_r2", fixed = TRUE)
})

test_that("the generated code compiles under Stan 2.32 (rstan)", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  for (label in names(specs)) {
    status <- rstan::stanc(model_code = .stan_code(specs[[label]]),
                           verbose = FALSE)$status
    expect_true(status, label = label)
  }
})

test_that("the generated code passes the syntax check of stanli", {
  skip_if_not_installed("stanli", "0.25.0")
  skip_if_not(stanli::stanli_available(), "the stanli runtime is not installed")
  for (label in names(specs)) {
    expect_true(stanli::stanli_check_syntax(.stan_code(specs[[label]])), label = label)
  }
})
