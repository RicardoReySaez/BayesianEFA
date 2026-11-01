# Bayesian fit indices and information criteria from the final draws

test_that("the sample moments count the means and the covariances", {
  expect_identical(.sample_moments("raw", 5), 20)
  expect_identical(.sample_moments("cov", 5), 15)
  expect_identical(.sample_moments("cor", 5), 10)
})

test_that("PSIS-LOO uses relative efficiencies for MCMC draws only", {
  # Autocorrelated draws, as MCMC gives, have relative efficiencies below 1.
  set.seed(8)
  ar <- function() as.numeric(stats::filter(stats::rnorm(100), 0.9, "recursive"))
  log_lik <- array(-2 + 0.1 * replicate(2 * 15, ar()), c(100, 2, 15))
  config <- validate_loo_args(list(r_eff = TRUE))
  with_eff <- .loo_fit(log_lik, "mcmc", config)
  without <- .loo_fit(log_lik, "mcmc", validate_loo_args(NULL))
  approx <- .loo_fit(log_lik, "advi", config)
  expect_s3_class(with_eff, "psis_loo")
  expect_equal(without$diagnostics, approx$diagnostics)
  expect_true(all(with_eff$diagnostics$n_eff < without$diagnostics$n_eff))
})

test_that("the indices follow the DevM formulas", {
  object <- .fake_befa()
  fm <- befa_fit_measures(object)
  expect_s3_class(fm, "befa_fitmeasures")
  expect_identical(posterior::variables(fm$fit_draws),
                   c("Chi2", "Chi2_Null", "BRMSEA", "BGamma", "Adj_BGamma",
                     "BMc", "SRMR", "BCFI", "BTLI"))
  expect_identical(unname(dim(fm$fit_draws)), c(100L, 2L, 9L))
  expect_named(fm$details, c("p_star", "pD", "N", "chi2_ppp"))

  d <- fm$posterior_fit
  N <- 60
  p_star <- 10
  df_H <- max(1, p_star - fm$details$pD)
  lambda_hat <- pmax(0, d$Chi2 - p_star)
  lambda_null <- pmax(0, d$Chi2_Null - p_star)
  expect_equal(d$BRMSEA, sqrt(lambda_hat / (df_H * N)))
  expect_equal(d$BGamma, 5 / (5 + 2 * lambda_hat / N))
  expect_equal(d$Adj_BGamma, 1 - (p_star / df_H) * (1 - d$BGamma))
  expect_equal(d$BMc, exp(-lambda_hat / (2 * N)))
  expect_equal(d$BCFI, pmin(1, pmax(0, 1 - lambda_hat / lambda_null)))
  expect_true(all(d$BTLI >= 0 & d$BTLI <= 1))
  expect_true(all(d$SRMR > 0 & d$SRMR < 0.3))
  expect_gte(fm$details$chi2_ppp, 0)
  expect_lte(fm$details$chi2_ppp, 1)
  expect_identical(fm$details$N, 60L)
})

test_that("the covariance and raw models compare with a fitted null model", {
  for (model in c("cov", "raw")) {
    object <- .fake_befa(model, missing = "FIML")
    null <- .null_draws(5, raw = model == "raw", n_iter = 100L)
    local_mocked_bindings(.run_stan = function(...) list(draws = null))
    fm <- befa_fit_measures(object, loo_args = list(r_eff = TRUE))
    expect_true(all(is.finite(as.matrix(fm$posterior_fit))))
    expect_identical(fm$details$p_star, .sample_moments(model, 5))
    expect_gt(min(fm$posterior_fit$Chi2_Null), max(fm$posterior_fit$Chi2))
  }
})

test_that("summary statistics and stored indices are handled", {
  object <- .fake_befa()
  object$fit_indices <- "stored"
  expect_identical(befa_fit_measures(object), "stored")
  object$fit_indices <- NULL
  object$data$summary <- TRUE
  expect_warning(res <- befa_fit_measures(object), "need the raw `data`")
  expect_null(res)
})

test_that("the fit measures and information criteria print as two tables", {
  fm <- befa_fit_measures(.fake_befa())
  out <- capture.output(res <- print(fm, digits = 2))
  expect_identical(res, fm)
  out <- paste(out, collapse = "\n")
  for (index in c("Table 1. Bayesian Fit Measures", "Chi2_ppp", "BRMSEA", "BTLI",
                  "Table 2. Information Criteria", "ELPD", "LOOIC", "DIC_i",
                  "BPIC", "SEs across persons")) {
    expect_match(out, index, fixed = TRUE)
  }
  tab <- build_fit_summary_table(fm)
  expect_identical(rownames(tab)[1:2], c("Chi2", "Chi2_ppp"))
  expect_identical(colnames(tab), c("Estimate", "SD", "CI_Low", "CI_High"))
  expect_true(all(is.na(tab["Chi2_ppp", -1])))
  expect_false("LOOIC" %in% rownames(tab))

  ic <- build_information_table(fm, probs = c(0.05, 0.95))
  looic <- fm$loo_object$estimates["looic", ]
  expect_equal(unname(ic["LOOIC", 1:2]), unname(looic))
  expect_equal(unname(ic["LOOIC", 3:4]),
               unname(looic[1] + stats::qnorm(c(0.05, 0.95)) * looic[2]))
})

test_that("the information criteria follow their definitions", {
  object <- .fake_befa("raw", missing = "FIML")
  m <- .posterior_metrics(object)
  lo <- .loo_fit(m$log_lik, "mcmc", validate_loo_args(NULL))
  ic <- .information_criteria(object, m$log_lik, lo)
  expect_identical(rownames(ic), c("ELPD", "LOOIC", "p_loo", "WAIC", "p_WAIC",
                                   "DIC", "p_D", "DIC_i", "p_V", "BPIC"))
  expect_true(all(ic[, "SE"] > 0))

  ll <- matrix(m$log_lik, prod(dim(m$log_lik)[1:2]))
  D <- -2 * rowSums(ll)
  # Deviance at the posterior means of the parameters, row by row
  p <- .fitted_parameters(object)
  Sigma <- tcrossprod(matrix(colMeans(p$lambda), 5, 2)) + diag(colMeans(p$psi))
  nu <- colMeans(p$nu)
  Y <- object$data$Y
  D_hat <- -2 * sum(vapply(seq_len(nrow(Y)), function(i) {
    o <- which(!is.na(Y[i, ]))
    mvnfast::dmvn(Y[i, o], nu[o], Sigma[o, o, drop = FALSE], log = TRUE)
  }, 0))
  p_D <- mean(D) - D_hat
  expect_equal(ic["p_D", "Estimate"], p_D)
  expect_equal(ic["DIC", "Estimate"], mean(D) + p_D)
  expect_equal(ic["BPIC", "Estimate"], mean(D) + 2 * p_D)
  expect_equal(ic["p_V", "Estimate"], stats::var(D) / 2)
  expect_equal(ic["DIC_i", "Estimate"], mean(D) + stats::var(D) / 2)
  waic <- suppressWarnings(loo::waic(ll))$estimates
  expect_equal(unname(ic["WAIC", ]), unname(waic["waic", ]))

  # The jackknife leaves each person out by recomputing the deviances.
  N <- ncol(ll)
  jack <- vapply(seq_len(N), function(i) {
    Di <- -2 * rowSums(ll[, -i, drop = FALSE])
    mean(Di) + stats::var(Di) / 2
  }, 0)
  expect_equal(ic["DIC_i", "SE"], sqrt((N - 1) / N * sum((jack - mean(jack))^2)))
})

test_that("the null model of an ADVI fit starts from one list of values", {
  object <- .fake_befa("cov")
  object$fit$settings <- .sampler_settings("advi", list(seed = 3))
  null <- .null_draws(5, raw = FALSE, n_iter = 200L)
  init <- NULL
  local_mocked_bindings(.run_stan = function(code, data, inits, ...) {
    init <<- inits
    list(draws = null)
  })
  res <- .null_metrics(object, ll_saturated = -400, loo_config = list(cores = 1))
  expect_named(init, "sigma")
  expect_length(res$chisq_vec, 400L)
})

test_that("a person without observed items adds nothing to the likelihood", {
  # Such rows are dropped when the data are prepared; the C++ code still
  # gives them a log-likelihood of zero.
  set.seed(10)
  Y <- rbind(matrix(stats::rnorm(60), 30, 2), c(NA, NA))
  res <- compute_posterior_metrics_cpp(
    lambda = matrix(c(0.6, 0.5), 1), psi = matrix(c(0.6, 0.7), 1),
    nu = matrix(0, 1, 2), M = 1L, n_iter = 1L, data = Y, M_sat = c(0, 0),
    S_sat = diag(2), ll_saturated = -80, model_type = 2L, has_missing = TRUE,
    Y_comp = Y[1:30, ], Y_miss = Y[31, , drop = FALSE], comp_idx_0 = 0:29,
    miss_idx_0 = 30L, obs_patterns_miss = list(integer(0)))
  expect_identical(res$log_lik[1, 1, 31], 0)
  Sigma <- tcrossprod(c(0.6, 0.5)) + diag(c(0.6, 0.7))
  expect_equal(res$log_lik[1, 1, 1:30],
               mvnfast::dmvn(Y[1:30, ], c(0, 0), Sigma, log = TRUE))
})
