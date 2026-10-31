test_that("the prior and the number of factors choose the loadings", {
  expect_identical(.befa_spec("cor", "unit_vector", 3L)$loadings, "uv")
  expect_identical(.befa_spec("cov", "unit_vector", 1L)$loadings, "uni")
  expect_identical(.befa_spec("raw", "normal", 1L)$loadings, "normal")
  spec <- .befa_spec("raw", "normal", 2L, patterns = "multiple")
  expect_s3_class(spec, "befa_spec")
  expect_identical(.spec_label(spec), "raw-normal-multiple-cholesky")
})

test_that("the normal prior is refused in the correlation model", {
  expect_error(.befa_spec("cor", "normal", 2L), "not available for the correlation")
})

test_that("the number of factors must be a positive whole number", {
  expect_error(.befa_spec("cor", "unit_vector", 0), "positive whole number")
  expect_error(.befa_spec("cor", "unit_vector", 1.5), "positive whole number")
  expect_error(.befa_spec("cor", "unit_vector", c(1, 2)), "positive whole number")
  expect_error(.befa_spec("cor", "unit_vector", NA), "positive whole number")
})

test_that("null models exist for the covariance and raw models only", {
  spec <- .befa_null_spec("raw", "multiple")
  expect_identical(.spec_label(spec), "raw-none-multiple-diagonal")
  expect_error(.befa_null_spec("cor"), "should be one of")
})

test_that("item SDs scale the unit-vector loadings outside the correlation model", {
  expect_false(.spec_has_scale(.befa_spec("cor", "unit_vector", 2L)))
  expect_true(.spec_has_scale(.befa_spec("cov", "unit_vector", 1L)))
  expect_false(.spec_has_scale(.befa_spec("raw", "normal", 2L)))
  expect_true(.spec_has_scale(.befa_null_spec("cov")))
})

test_that("matrices are checked", {
  expect_error(.check_matrix(1:3, "x"), "`x` must be a numeric matrix")
  expect_error(.check_matrix(matrix("a"), "x"), "numeric matrix")
  expect_true(.check_matrix(diag(2), "x"))
})
