# The C++ factor scores and posterior metrics are still the v0.1.0 code. These
# tests pin their output, random draws included, before they are adapted to
# the new parameterization.

for (case in c("factor_scores_cor_complete", "factor_scores_cor_missing")) {
  test_that(paste("compute_factor_scores_cpp() reproduces", case), {
    fx <- .read_fixture(case)
    set.seed(fx$inputs$seed)
    args <- fx$inputs[setdiff(names(fx$inputs), "seed")]
    expect_equal(do.call(compute_factor_scores_cpp, args), fx$output,
                 tolerance = 1e-10)
  })
}

metric_cases <- c(posterior_metrics_cor_complete = "cor",
                  posterior_metrics_raw_complete = "raw",
                  posterior_metrics_raw_missing = "raw")

for (case in names(metric_cases)) {
  test_that(paste("compute_posterior_metrics_cpp() reproduces", case), {
    fx <- .read_fixture(case)
    set.seed(fx$inputs$seed)
    args <- fx$inputs[setdiff(names(fx$inputs), "seed")]
    expect_equal(do.call(compute_posterior_metrics_cpp, args), fx$output,
                 tolerance = 1e-10)
  })

  test_that(paste("compute_posterior_metrics() reproduces", case), {
    fx <- .read_fixture(case)
    model_type <- metric_cases[[case]]
    J <- ncol(fx$inputs$data)
    cells <- paste0("[", rep(seq_len(J), J), ",", rep(seq_len(J), each = J), "]")
    vars <- c(paste0(if (model_type == "cor") "Rho" else "Sigma", cells),
              if (model_type == "raw") paste0("nu[", seq_len(J), "]"))
    draws <- fx$inputs$post_draws
    dimnames(draws) <- list(NULL, NULL, vars)
    # The R wrapper only reads the draws and the model type, and
    # as_draws_array() leaves a draws_array as it is.
    object <- list(model_type = model_type,
                   stanfit = posterior::as_draws_array(draws))

    set.seed(fx$inputs$seed)
    res <- compute_posterior_metrics(object, fx$inputs$data)
    expect_equal(res[names(fx$output)], fx$output, tolerance = 1e-10)
    expect_equal(res$ll_saturated, fx$inputs$ll_saturated, tolerance = 1e-10)
  })
}
