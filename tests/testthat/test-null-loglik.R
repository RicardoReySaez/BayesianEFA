test_that("the null log-likelihood matches v0.1.0 FIML with a diagonal covariance", {
  fx <- .read_fixture("fiml_loglik_null_diag")
  Y <- .fiml_data_matrix(fx$inputs$data_fiml)
  ll <- null_loglik_cpp(Y, matrix(fx$inputs$mu, nrow = 1),
                        matrix(sqrt(diag(fx$inputs$Sigma)), nrow = 1))
  expect_equal(as.vector(ll), fx$output, tolerance = 1e-10)
})

test_that("a person with no observed item adds nothing", {
  Y <- rbind(c(0.5, -1), c(NA, NA))
  ll <- null_loglik_cpp(Y, matrix(0, 2, 2), matrix(1, 2, 2))
  expect_identical(ll[, 2], c(0, 0))
  expect_equal(ll[, 1], rep(sum(stats::dnorm(c(0.5, -1), log = TRUE)), 2))
})

test_that("draws of the wrong shape are refused", {
  Y <- matrix(0, 3, 2)
  expect_error(null_loglik_cpp(Y, matrix(0, 4, 2), matrix(1, 4, 3)), "S x J")
  expect_error(null_loglik_cpp(Y, matrix(0, 5, 2), matrix(1, 4, 2)), "S x J")
})

# This function returns fixed null-model draws in place of a Stan fit: two
# chains of five draws of sigma, and of nu for the raw model.
.null_draws <- function(J, raw, seed = 21) {
  set.seed(seed)
  vars <- c(paste0("sigma[", seq_len(J), "]"),
            if (raw) paste0("nu[", seq_len(J), "]"))
  values <- c(stats::runif(5 * 2 * J, 0.8, 1.6),
              if (raw) stats::rnorm(5 * 2 * J, 0, 0.3))
  posterior::as_draws_array(array(values, dim = c(5, 2, length(vars)),
                                  dimnames = list(NULL, NULL, vars)))
}

for (model_type in c("raw", "cov")) {
  test_that(paste("compute_null_metrics() keeps the v0.1.0 chi-square for the",
                  model_type, "model"), {
    set.seed(22)
    J <- 4
    Y <- matrix(stats::rnorm(30 * J, 1, 1.2), ncol = J)
    Y[c(2, 7, 15, 23)] <- NA
    draws <- .null_draws(J, raw = model_type == "raw")
    local_mocked_bindings(fit_befa_model = function(...) draws)

    # The raw case also exercises the relative efficiencies of loo.
    res <- compute_null_metrics(
      data = Y, model_type = model_type, ll_saturated = -150,
      stan_data = list(), stan_args_list = list(), lambda_prior = "unit_vector",
      loo_config = list(r_eff = model_type == "raw", cores = 1)
    )

    # Reference: the R loop of v0.1.0 over draws, chain after chain
    data <- if (model_type == "cov") scale(Y, center = TRUE, scale = FALSE) else Y
    data_fiml <- list(
      comp_idx = which(stats::complete.cases(data)),
      miss_idx = which(!stats::complete.cases(data)),
      Y_comp = data[stats::complete.cases(data), , drop = FALSE],
      Y_miss = data[!stats::complete.cases(data), , drop = FALSE]
    )
    expected <- numeric(0)
    for (ch in 1:2) for (i in 1:5) {
      sigma <- as.numeric(draws[i, ch, paste0("sigma[", 1:J, "]")])
      nu <- if (model_type == "raw") {
        as.numeric(draws[i, ch, paste0("nu[", 1:J, "]")])
      } else {
        rep(0, J)
      }
      ll <- compute_fiml_loglik(data_fiml, mu = nu, Sigma = diag(sigma^2))
      expected <- c(expected, 2 * (-150 - sum(ll)))
    }

    p_star <- if (model_type == "raw") J * (J + 3) / 2 else J * (J + 1) / 2
    expect_equal(res$chisq_vec, expected, tolerance = 1e-10)
    expect_equal(res$lambda_vec, pmax(0, expected - p_star), tolerance = 1e-10)
    expect_gte(res$df, 1)
  })
}
