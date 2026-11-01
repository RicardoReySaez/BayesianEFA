cases <- c("rsp_align_J9_M3", "rsp_align_J6_M2", "rsp_align_J5_M1")
formats <- c("column_major", "row_major")

for (case in cases) for (fmt in formats) {
  test_that(paste("rsp_align() reproduces v0.1.0 for", case, "in", fmt, "format"), {
    fx <- .read_fixture(paste0(case, "_", fmt))
    res <- do.call(rsp_align, fx$inputs)
    res$summary <- NULL
    # The per-draw Varimax now runs in C++; the rest of the algorithm is the
    # published one, so only rounding may differ. The algorithm stops when the
    # objective stops improving, so once it is flat a rounding-size rise can
    # add one iteration on some platforms: the trace is compared by its final
    # value.
    expect_equal(res[names(res) != "objective"],
                 fx$output[names(fx$output) != "objective"], tolerance = 1e-8)
    expect_equal(tail(res$objective, 1L), tail(fx$output$objective, 1L),
                 tolerance = 1e-8)
    expect_lte(abs(length(res$objective) - length(fx$output$objective)), 1L)
  })
}

test_that("rsp_align() accepts posterior draws_matrix objects", {
  fx <- .read_fixture("rsp_align_J6_M2_column_major")
  draws <- posterior::as_draws_matrix(fx$inputs$lambda_draws)
  res <- rsp_align(draws, n_items = 6, n_factors = 2, n_chains = 2,
                   format = "column_major")
  expect_equal(unname(res$Lambda_hat_mcmc), unname(fx$output$Lambda_hat_mcmc),
               tolerance = 1e-8)
})
