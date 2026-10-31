# The C++ factor scores and posterior metrics take the loadings,
# uniquenesses and means of each draw (chains stacked) instead of the v0.1.0
# draws of Sigma. These tests feed them the v0.1.0 fixtures, translated to
# the new inputs, and expect the same output, random draws included.

# This function returns some variables (0-based) of a v0.1.0 draws array
# (iterations x chains x variables) with the chains stacked.
.stack_fixture <- function(draws, idx) {
  matrix(draws[, , idx + 1L], nrow = prod(dim(draws)[1:2]))
}

# This function returns, draw by draw, a lower-triangular root of each
# covariance matrix in a stacked matrix of them (one row per draw), as
# loadings with as many factors as items.
.root_draws <- function(sigma, J) {
  t(apply(sigma, 1L, function(s) as.vector(t(chol(matrix(s, J, J))))))
}

for (case in c("factor_scores_cor_complete", "factor_scores_cor_missing")) {
  test_that(paste("compute_factor_scores_cpp() reproduces", case), {
    fx <- .read_fixture(case)
    inp <- fx$inputs
    J <- inp$J
    M <- inp$M
    lambda <- .stack_fixture(inp$post_draws, inp$lambda_idx)
    sigma <- .stack_fixture(inp$post_draws, inp$sig_idx)
    psi <- sigma[, (seq_len(J) - 1L) * J + seq_len(J)] - .item_sums(lambda^2, J, M)
    set.seed(inp$seed)
    eta <- compute_factor_scores_cpp(lambda, psi, matrix(0, nrow(psi), J),
                                     inp$data, J, M, inp$has_missing,
                                     inp$obs_patterns)
    expect_equal(eta, matrix(fx$output, nrow = nrow(lambda)), tolerance = 1e-10)
  })
}

metric_cases <- c(posterior_metrics_cor_complete = "cor",
                  posterior_metrics_raw_complete = "raw",
                  posterior_metrics_raw_missing = "raw")

for (case in names(metric_cases)) {
  test_that(paste("compute_posterior_metrics_cpp() reproduces", case), {
    fx <- .read_fixture(case)
    inp <- fx$inputs
    J <- ncol(inp$data)
    S <- prod(dim(inp$post_draws)[1:2])
    lambda <- .root_draws(.stack_fixture(inp$post_draws, inp$sig_idx), J)
    nu <- if (length(inp$nu_idx)) {
      .stack_fixture(inp$post_draws, inp$nu_idx)
    } else {
      matrix(0, S, J)
    }
    args <- inp[setdiff(names(inp), c("seed", "post_draws", "sig_idx", "nu_idx"))]
    args <- c(list(lambda = lambda, psi = matrix(0, S, J), nu = nu, M = J,
                   n_iter = dim(inp$post_draws)[1L]), args)
    set.seed(inp$seed)
    expect_equal(do.call(compute_posterior_metrics_cpp, args), fx$output,
                 tolerance = 1e-10)
  })

  test_that(paste(".posterior_metrics() reproduces", case), {
    fx <- .read_fixture(case)
    inp <- fx$inputs
    model_type <- metric_cases[[case]]
    J <- ncol(inp$data)
    S <- prod(dim(inp$post_draws)[1:2])
    n_iter <- dim(inp$post_draws)[1L]
    lambda <- .root_draws(.stack_fixture(inp$post_draws, inp$sig_idx), J)
    std <- model_type == "cor"
    parts <- list(lambda, matrix(0, S, J))
    names(parts) <- if (std) c("Lambda_std", "Psi_std") else c("Lambda", "Psi")
    if (length(inp$nu_idx)) parts$nu <- .stack_fixture(inp$post_draws, inp$nu_idx)
    # The R wrapper reads the data, the model type, the dimensions and the
    # final draws of the fit.
    object <- list(model_type = model_type, data = list(Y = inp$data),
                   stan_data = list(J = J), n_factors = J,
                   draws = .as_chains(.named_draws(parts, J, J), n_iter, S / n_iter))

    set.seed(inp$seed)
    res <- .posterior_metrics(object)
    expect_equal(res[names(fx$output)], fx$output, tolerance = 1e-10)
    expect_equal(res$ll_saturated, inp$ll_saturated, tolerance = 1e-10)
  })
}
