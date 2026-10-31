test_that("errors and warnings join their pieces and hide the call", {
  err <- tryCatch(.befa_stop("`n_factors` must be ", 2L, "."),
                  error = function(e) e)
  expect_identical(conditionMessage(err), "`n_factors` must be 2.")
  expect_null(conditionCall(err))

  wrn <- tryCatch(.befa_warn("three ", "pieces"), warning = function(w) w)
  expect_identical(conditionMessage(wrn), "three pieces")
  expect_null(conditionCall(wrn))
})

test_that("messages name the package", {
  expect_message(.befa_msg("using ", "stanli."), "^BayesianEFA: using stanli\\.")
})
