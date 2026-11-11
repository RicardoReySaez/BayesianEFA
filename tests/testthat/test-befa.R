# Argument checks: they all stop before any sampling.

test_that("befa() checks where the data come from", {
  Y <- .simulated_data()
  expect_error(befa(Y, 2, sample_nobs = 40, sample_cov = cov(Y)), "not both")
  expect_error(befa(n_factors = 2), "Provide `data`, or summary statistics")
  expect_error(befa(sample_nobs = 40, sample_cov = cov(Y), n_factors = 2,
                    missing = "FIML"), "FIML needs the raw `data`")
  expect_error(befa(list(a = 1:3), 1), "data frame or a matrix")
  expect_error(befa(data.frame(a = letters[1:5], b = 1:5), 1), "numeric columns")
  expect_error(befa(sample_nobs = 40, sample_cov = cov(Y), sample_cor = cor(Y),
                    n_factors = 2), "either `sample_cov` or `sample_cor`")
  expect_error(befa(sample_nobs = 40, sample_cov = cov(Y)[, -1], n_factors = 2),
               "symmetric square matrix")
})

test_that("befa() checks the model it is asked to fit", {
  Y <- .simulated_data()
  expect_error(befa(Y), "`n_factors` is required")
  expect_error(befa(Y, 0), "positive whole number")
  expect_error(befa(Y, 5), "smaller than the number of items \\(5\\)")
  expect_error(befa(Y, 2, ordered = TRUE), "not implemented yet")
  expect_error(befa(Y, 2, verbose = "yes"), "`verbose` must be TRUE or FALSE")
  expect_error(befa(Y, 2, model = "cor", lambda_prior = "normal"),
               "not available for the correlation model")
  expect_error(befa(Y, 2, chainz = 2), "`chainz` is not an option of MCMC")
  expect_error(befa(Y, 2, backend = "brms"), "must be \"auto\" or one of")
  expect_error(befa(Y, 2, prior_only = NA),
               "`prior_only` must be TRUE or FALSE")
  # Asking for factor scores without the likelihood warns before any sampling
  # (the number of factors then stops the fit).
  expect_warning(
    expect_error(befa(Y, 5, prior_only = TRUE, factor_scores = TRUE),
                 "smaller than the number of items"),
    "Factor scores need the likelihood")
})

test_that("starting values are built, reused or checked", {
  spec <- .befa_spec("cor", "unit_vector", 2L)
  prep <- .prepared_for(spec)
  prior <- .prior_for(spec)
  mcmc <- .sampler_settings("mcmc", list(chains = 3, seed = 1))
  expect_length(.befa_inits(mcmc, spec, prep, 2L, prior), 3L)
  one <- .start_values(spec, prep, 2L, prior, 1L, 1L)[[1L]]
  mcmc$init <- one
  expect_identical(.befa_inits(mcmc, spec, prep, 2L, prior), rep(list(one), 3))
  mcmc$init <- list(one, one)
  expect_error(.befa_inits(mcmc, spec, prep, 2L, prior), "has 2 lists")
  mcmc$init <- 0.5
  expect_error(.befa_inits(mcmc, spec, prep, 2L, prior), "must be a list")
  advi <- .sampler_settings("advi", list(seed = 1))
  expect_named(.befa_inits(advi, spec, prep, 2L, prior), c("h2", "Z"))
  path <- .sampler_settings("pathfinder", list(seed = 1, num_paths = 2))
  expect_length(.befa_inits(path, spec, prep, 2L, prior), 2L)
})

# Fits with stanli, which needs no compilation. They are short, so the
# convergence warnings they raise are expected and silenced.

.skip_without_stanli <- function() {
  skip_on_cran()
  skip_if_not_installed("stanli", "0.25.0")
  skip_if_not(stanli::stanli_available(), "the stanli runtime is not installed")
}

# This function fits a short stanli model quietly: two MCMC chains of 150
# draws, or the defaults of ADVI and Pathfinder.
.quick_befa <- function(data = .simulated_data(N = 300), n_factors = 2L,
                        method = "mcmc", ...) {
  withr::local_options(BayesianEFA.cache_dir = withr::local_tempdir())
  short <- if (method == "mcmc") {
    list(chains = 2, iter_warmup = 200, iter_sampling = 150)
  }
  args <- c(list(data, n_factors, method = method, backend = "stanli",
                 seed = 5, verbose = FALSE), short, list(...))
  suppressWarnings(do.call(befa, args))
}

# This function returns the posterior means of the standardized loadings.
.mean_loadings <- function(fit, M = fit$n_factors) {
  draws <- extract_posterior_draws(fit, "Lambda_std", metric = "std")
  matrix(colMeans(draws), ncol = M)
}

# The true loadings of .simulated_data(): items 1-3 on one factor, 4-5 on
# the other.
.true_pattern <- matrix(c(1, 1, 1, 0, 0, 0, 0, 0, 1, 1), 5, 2)

grid <- expand.grid(model = c("cor", "cov", "raw"),
                    lambda_prior = c("unit_vector", "normal"),
                    missing = c("listwise", "FIML"), stringsAsFactors = FALSE)
grid <- grid[!(grid$model == "cor" & grid$lambda_prior == "normal"), ]

for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  label <- paste(g$model, g$lambda_prior, g$missing, sep = " / ")
  test_that(paste("befa() fits and aligns", label), {
    .skip_without_stanli()
    data <- .simulated_data(N = 300, missing = g$missing == "FIML")
    fit <- .quick_befa(data, model = g$model, lambda_prior = g$lambda_prior,
                       missing = g$missing)
    expect_s3_class(fit, "befa")
    expect_identical(unname(dim(fit$draws)[1:2]), c(150L, 2L))
    expect_identical(fit$spec$patterns, if (g$missing == "FIML") "multiple" else "single")

    vars <- unique(sub("\\[.*", "", posterior::variables(fit$draws)))
    expected <- if (g$model == "cor") {
      c("Lambda_std", "Psi_std", "h2")
    } else if (g$lambda_prior == "normal") {
      c("Lambda", "Psi", if (g$model == "raw") "nu")
    } else {
      c("Lambda", "Psi", "sigma", "h2", if (g$model == "raw") "nu")
    }
    expect_identical(vars, expected)

    # Each true factor is recovered by one aligned column with positive mean.
    L <- .mean_loadings(fit)
    hit <- abs(stats::cor(L, .true_pattern))
    expect_gt(min(apply(hit, 2, max)), 0.9)
    expect_true(all(colMeans(L) > 0))
    # F1 explains the most variance.
    expect_gte(sum(L[, 1]^2), sum(L[, 2]^2))
    std <- extract_posterior_draws(fit, c("Psi_std", "h2"), metric = "std")
    expect_true(all(std > 0 & std < 1))
    expect_named(fit$rsp, c("Lambda_star", "objective", "sign_vectors",
                            "perm_vectors", "factor_order"))
    expect_false(is.na(fit$diagnostics$max_rhat))
  })
}

test_that("befa() fits one factor and summary statistics", {
  .skip_without_stanli()
  Y <- .simulated_data(N = 300)
  uni <- .quick_befa(Y[, 1:3], n_factors = 1L)
  expect_identical(posterior::variables(uni$draws)[1:3],
                   c("Lambda_std[1,1]", "Lambda_std[2,1]", "Lambda_std[3,1]"))
  expect_true(all(.mean_loadings(uni) > 0.4))

  withr::local_options(BayesianEFA.cache_dir = withr::local_tempdir())
  summ <- suppressWarnings(befa(sample_nobs = 300, sample_cov = stats::cov(Y),
                                sample_mean = colMeans(Y), n_factors = 2,
                                model = "raw", backend = "stanli", chains = 2,
                                iter_warmup = 200, iter_sampling = 150, seed = 5,
                                verbose = FALSE))
  expect_true(summ$data$summary)
  expect_equal(unname(colMeans(extract_posterior_draws(summ, "nu"))), colMeans(Y),
               tolerance = 0.1)
})

test_that("both likelihood engines give the same posterior", {
  .skip_without_stanli()
  data <- .simulated_data(N = 300, missing = TRUE)
  woodbury <- .quick_befa(data, model = "raw", missing = "FIML",
                          compute_fit_indices = FALSE)
  withr::local_options(BayesianEFA.engine = "cholesky")
  cholesky <- .quick_befa(data, model = "raw", missing = "FIML",
                          compute_fit_indices = FALSE)
  expect_identical(woodbury$spec$engine, "woodbury")
  expect_identical(cholesky$spec$engine, "cholesky")
  mean_of <- function(fit, par) colMeans(.stacked_draws(fit$draws, par))
  expect_lt(max(abs(mean_of(woodbury, "h2") - mean_of(cholesky, "h2"))), 0.05)
  expect_lt(max(abs(mean_of(woodbury, "nu") - mean_of(cholesky, "nu"))), 0.1)
})

test_that("befa() runs ADVI and Pathfinder and can skip the alignment", {
  .skip_without_stanli()
  advi <- .quick_befa(method = "advi")
  expect_identical(advi$method, "advi")
  expect_true("lp_approx__" %in% posterior::variables(advi$fit$draws))
  path <- .quick_befa(method = "pathfinder", rotate = "none")
  expect_null(path$rsp)
  expect_identical(path$rotation, "none")
})

test_that("befa() reports the fit and warns when it is too short", {
  .skip_without_stanli()
  withr::local_options(BayesianEFA.cache_dir = withr::local_tempdir())
  msgs <- character(0)
  expect_warning(
    withCallingHandlers(
      befa(.simulated_data(N = 300), 2, backend = "stanli", chains = 2,
           iter_warmup = 100, iter_sampling = 50, seed = 3),
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
      }),
    "may be unreliable")
  msgs <- paste(msgs, collapse = "")
  expect_match(msgs, "Bayesian exploratory factor analysis")
  expect_match(msgs, "Sampled 2 chains")
  expect_match(msgs, "Largest R-hat")
  expect_match(msgs, "Computing the fit indices")
})

test_that("befa() adds the scores, omega and fit indices of a fit", {
  .skip_without_stanli()
  fit <- .quick_befa(model = "cov", factor_scores = TRUE)
  expect_s3_class(fit$eta, "draws_array")
  expect_identical(posterior::nvariables(fit$eta), 600L)
  expect_s3_class(fit$reliability, "befa_reliability")
  expect_s3_class(fit$fit_indices, "befa_fitmeasures")
  # The null covariance model is fitted with stanli and fits far worse.
  fm <- fit$fit_indices$posterior_fit
  expect_gt(min(fm$Chi2_Null), max(fm$Chi2))
  expect_gt(mean(fm$BCFI), 0.9)
  expect_true("eta[1,1]" %in% posterior::variables(
    extract_posterior_draws(fit, format = "array", metric = "std")))
  out <- utils::capture.output(print(summary(fit)))
  expect_true(any(grepl("Table 4", out)))

  again <- .quick_befa(model = "cov", factor_scores = TRUE,
                       compute_fit_indices = FALSE, compute_reliability = FALSE)
  expect_identical(again$eta, fit$eta)
  expect_null(again$fit_indices)
  expect_null(again$reliability)
})

test_that("the scores and fit indices are skipped with summary statistics", {
  object <- .fake_befa()
  object$data$summary <- TRUE
  expect_warning(
    expect_warning(res <- .befa_downstream(object, TRUE, TRUE, TRUE, 1L, FALSE),
                   "Factor scores need the raw `data`"),
    "Fit indices need the raw `data`")
  expect_null(res$eta)
  expect_null(res$fit_indices)
  expect_s3_class(res$reliability, "befa_reliability")
})

test_that("the scores follow the seed and a failing fit index only warns", {
  object <- .fake_befa()
  expect_message(first <- .befa_downstream(object, TRUE, FALSE, FALSE, 11L, TRUE),
                 "factor scores")
  second <- .befa_downstream(object, TRUE, FALSE, FALSE, 11L, FALSE)
  expect_identical(second$eta, first$eta)
  expect_null(first$reliability)

  local_mocked_bindings(.fit_measures = function(...) stop("singular"))
  expect_warning(res <- .befa_downstream(object, FALSE, FALSE, TRUE, 1L, FALSE),
                 "could not be computed: singular")
  expect_null(res$fit_indices)
})

test_that("befa() draws from the prior alone with prior_only = TRUE", {
  .skip_without_stanli()
  withr::local_options(BayesianEFA.cache_dir = withr::local_tempdir())
  msgs <- character(0)
  fit <- suppressWarnings(withCallingHandlers(
    befa(.simulated_data(N = 300), 2, prior_only = TRUE, backend = "stanli",
         chains = 2, iter_warmup = 200, iter_sampling = 150, seed = 5,
         refresh = 0),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }))
  expect_match(paste(msgs, collapse = ""), "Sampling from the prior alone")
  expect_true(fit$options$prior_only)
  expect_identical(fit$stan_data$prior_only, 1L)
  expect_null(fit$fit_indices)
  expect_s3_class(fit$reliability, "befa_reliability")
  # The default Beta(1, 1) prior makes each communality uniform on (0, 1),
  # with mean 0.5 and SD 0.29, while their posterior from 300 rows has an SD
  # near 0.05. With 300 draws of 5 items, the mean is well within 0.1 of 0.5.
  h2 <- as.vector(.stacked_draws(fit$draws, "h2"))
  expect_lt(abs(mean(h2) - 0.5), 0.1)
  expect_gt(stats::sd(h2), 0.2)

  fits <- .quick_befa(n_factors = 1:2, prior_only = TRUE)
  expect_s3_class(fits, "befa_list")
  expect_null(attr(fits, "comparison"))
  expect_true(all(vapply(fits, .is_prior_only, TRUE)))
  expect_output(print(fits), "prior alone")
})

test_that("what needs the likelihood refuses a fit from the prior alone", {
  object <- .fake_befa()
  object$options$prior_only <- TRUE
  expect_warning(expect_null(befa_fit_measures(object)), "prior alone")
  expect_error(befa_factor_scores(object), "Factor scores need the likelihood")
  expect_error(befa_bridge_sampler(object), "needs draws from the posterior")
  expect_error(befa_compare(object, .fake_befa()),
               "Model 1 drew from the prior")
  out <- utils::capture.output(print(summary(object)))
  expect_true(any(grepl("Draws from the prior alone", out)))
  expect_true(any(grepl("Prior means", out)))
  # Fits stored before `prior_only` existed drew from the posterior.
  object$options$prior_only <- NULL
  expect_false(.is_prior_only(object))
})

test_that("befa() reproduces v0.1.0 on the Holzinger-Swineford data", {
  .skip_without_stanli()
  withr::local_options(BayesianEFA.cache_dir = withr::local_tempdir())
  fit <- befa(HS_data, n_factors = 3, backend = "stanli", seed = 17,
              verbose = FALSE)
  # Posterior means printed in the v0.1.0 documentation of befa().
  v010_h2 <- c(0.476, 0.240, 0.447, 0.715, 0.751, 0.687, 0.486, 0.527, 0.442)
  v010_L <- matrix(c(0.594, 0.459, 0.648, 0.110, 0.029, 0.156, -0.053, 0.171, 0.401,
                     0.313, 0.130, 0.076, 0.833, 0.862, 0.809, 0.099, 0.075, 0.165,
                     0.122, -0.038, 0.108, 0.073, 0.068, 0.062, 0.681, 0.693, 0.492),
                   9, 3)
  expect_equal(unname(colMeans(extract_posterior_draws(fit, "h2"))), v010_h2,
               tolerance = 0.02)
  # The order of the factors is arbitrary, so each v0.1.0 column is matched
  # with the closest aligned column.
  L <- .mean_loadings(fit)
  for (m in 1:3) {
    best <- which.max(abs(stats::cor(v010_L[, m], L)))
    expect_lt(max(abs(L[, best] - v010_L[, m])), 0.05)
  }
  expect_length(fit$diagnostics$issues, 0L)
})

test_that("extract_posterior_draws() and posterior_summaries() read befa fits", {
  .skip_without_stanli()
  fit <- .quick_befa(model = "cov")
  expect_s3_class(extract_posterior_draws(fit, "Lambda"), "draws_matrix")
  expect_s3_class(extract_posterior_draws(fit, format = "rvars"), "draws_rvars")
  expect_s3_class(extract_posterior_draws(fit, format = "df"), "draws_df")
  expect_s3_class(extract_posterior_draws(fit, format = "list"), "draws_list")
  expect_s3_class(extract_posterior_draws(fit, format = "array"), "draws_array")
  summ <- posterior_summaries(fit, "h2", metric = "std")
  expect_identical(summ$variable, paste0("h2[", 1:5, "]"))
  expect_error(extract_posterior_draws(list()), "fit returned by befa")
})
