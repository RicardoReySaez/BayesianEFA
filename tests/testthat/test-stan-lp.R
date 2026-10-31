# These tests evaluate the log density of every generated program with stanli,
# which needs no compilation, and compare it with a reference computed in R
# row by row. Agreement checks the code, the scatter roots, the patterns and
# every additive constant at once.

test_that("the repulsion constant normalizes ||z||^xi exp(-||z||^2 / 2)", {
  for (M in 2:4) for (xi in c(0, 3, 30)) {
    surface <- 2 * pi^(M / 2) / gamma(M / 2)
    radial <- stats::integrate(function(r) r^(xi + M - 1) * exp(-r^2 / 2),
                               0, Inf, rel.tol = 1e-12)$value
    expect_equal(.repulsion_log_const(xi, M), log(surface * radial),
                 tolerance = 1e-10)
  }
  # Without repulsion it is the constant of the standard normal.
  expect_equal(.repulsion_log_const(0, 3), 1.5 * log(2 * pi))
})

# This function returns the log density that stanli gives for the program of
# a structure, at constrained values of its parameters.
.stanli_lp <- function(spec, data, pars) {
  m <- stanli::stanli_model(code = .stan_code(spec), data = data)
  stanli::log_prob_grad(m, stanli::unconstrain(m, pars))$lp
}

for (label in names(.all_specs())) {
  test_that(paste("the log density of", label, "matches the reference"), {
    skip_if_not_installed("stanli", "0.25.0")
    skip_if_not(stanli::stanli_available(), "the stanli runtime is not installed")
    spec <- .all_specs()[[label]]
    prep <- .prepared_for(spec)
    prior <- .prior_for(spec)
    M <- .spec_factors(spec)
    data <- .stan_data(spec, prep, M, prior)
    pars <- .random_pars(spec, prep$J, M)
    item_sd <- if (spec$loadings == "normal") .item_sd(prep)
    expect_equal(.stanli_lp(spec, data, pars),
                 .reference_lp(spec, prep, pars, prior, item_sd),
                 tolerance = 1e-8)
  })
}

test_that("summary statistics give the same log density as the rows", {
  skip_if_not_installed("stanli", "0.25.0")
  skip_if_not(stanli::stanli_available(), "the stanli runtime is not installed")
  Y <- .simulated_data()
  for (model_type in c("cor", "cov", "raw")) {
    spec <- .befa_spec(model_type, "unit_vector", 2L)
    prior <- .prior_for(spec)
    pars <- .random_pars(spec, ncol(Y), 2L)
    rows <- .prepare_befa_data(Y, model_type)
    summ <- .prepare_befa_data(model_type = model_type, sample_nobs = nrow(Y),
                               sample_cov = stats::cov(Y), sample_mean = colMeans(Y))
    lp_rows <- .stanli_lp(spec, .stan_data(spec, rows, 2L, prior), pars)
    lp_summ <- .stanli_lp(spec, .stan_data(spec, summ, 2L, prior), pars)
    expect_equal(lp_summ, lp_rows, tolerance = 1e-8, label = model_type)
    expect_equal(lp_summ, .reference_lp(spec, summ, pars, prior), tolerance = 1e-8)
  }
})
