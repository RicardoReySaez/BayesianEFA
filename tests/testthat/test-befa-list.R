# One model per number of factors: befa(n_factors = 1:3)

# This function returns a fake fit of M factors with its omega coefficients
# and fit indices, as befa() stores them.
.fake_fit <- function(M, ...) {
  object <- .fake_befa(M = M, ...)
  object$reliability <- befa_reliability(object)
  object$fit_indices <- befa_fit_measures(object)
  object
}

# This function replaces the Stan fit of each model by a fake one, in the
# calling test, and records the calls of each fit. The models in `fail` stop
# with an error.
.mock_fits <- function(fail = integer(0), env = parent.frame()) {
  calls <- new.env(parent = emptyenv())
  calls$run <- list()
  testthat::local_mocked_bindings(
    .resolve_backend = function(...) "stanli",
    .befa_fit_one = function(inputs, run, shared = NULL) {
      if (inputs$n_factors %in% fail) stop("did not converge")
      calls$run[[length(calls$run) + 1L]] <- run
      calls$shared <- c(calls$shared, list(shared))
      .fake_fit(inputs$n_factors)
    },
    .env = env
  )
  calls
}

# This function runs befa() on the data of .fake_befa() with fake fits.
.fake_list <- function(n_factors = 1:3, ...) {
  # Three factors for five items are above the Ledermann bound.
  withCallingHandlers(
    befa(.simulated_data(N = 60L, seed = 7), n_factors, verbose = FALSE, ...),
    warning = function(w) {
      if (grepl("Ledermann", conditionMessage(w))) invokeRestart("muffleWarning")
    })
}

# This function returns what print() shows as one line.
.shown <- function(x, ...) {
  gsub("\\s+", " ", paste(utils::capture.output(print(x, ...)), collapse = " "))
}

test_that("the numbers of factors are checked, sorted and made unique", {
  expect_identical(.factor_numbers(c(3, 1, 2, 1), 9L), 1:3)
  expect_identical(.factor_numbers(2, 9L), 2L)
  expect_error(.factor_numbers(c(1, 5), 5L), "smaller than the number of items \\(5\\)")

  # (J - M)^2 < J + M: with 6 items, 4 factors have more parameters than moments.
  expect_warning(.factor_numbers(c(3, 4), 6L), "With 4 factors and 6 items")
  expect_no_warning(.factor_numbers(c(1, 3), 6L))
  expect_no_warning(.factor_numbers(4, 6L))

  Y <- .simulated_data()
  expect_error(befa(Y, c(1, 1.5)), "positive whole numbers")
  expect_error(befa(Y, c(0, 2)), "positive whole numbers")
  expect_error(befa(Y, c(1, NA)), "positive whole numbers")
  expect_error(befa(Y, numeric(0)), "positive whole numbers")
})

test_that("the comparison arguments and the starting values are checked", {
  local_mocked_bindings(.resolve_backend = function(...) "stanli")
  Y <- .simulated_data()
  expect_error(befa(Y, 1:2, compare = "yes"), "`compare` must be TRUE or FALSE")
  expect_error(befa(Y, 1:2, compare_criteria = c("fit", "bic")), "`compare_criteria`")
  expect_error(befa(Y, 1:2, compare_criteria = 1), "`compare_criteria`")
  expect_error(befa(Y, 1:2, init = list(h2 = 0.5)), "`init` depends on")
})

test_that("each number of factors takes the priors that apply to it", {
  given <- list(xi = 50, lambda = c(2, 2))
  expect_no_warning(priors <- .befa_priors(given, "cor", "unit_vector", 1:2,
                                           5L))
  expect_identical(priors[[1L]]$lambda, c(2, 2))
  expect_null(priors[[1L]]$xi)
  expect_identical(priors[[2L]]$xi, 50)
  expect_null(priors[[2L]]$lambda)

  # A name that applies to no model is reported once.
  n <- 0L
  withCallingHandlers(.befa_priors(list(foo = 1), "cor", "unit_vector", 1:3,
                                   5L),
                      warning = function(w) {
                        n <<- n + 1L
                        expect_match(conditionMessage(w), "foo")
                        invokeRestart("muffleWarning")
                      })
  expect_identical(n, 1L)
  # One number of factors keeps the warning about the other kind of model.
  expect_warning(.befa_priors(list(xi = 50), "cor", "unit_vector", 1L, 5L),
                 "xi")
  expect_identical(.befa_priors(list(), "cor", "unit_vector", 2L, 5L)[[1L]],
                   prepare_befa_priors(list(), "cor", "unit_vector", 2L, 5L))
  # The default prior on the communalities follows each number of factors.
  priors <- .befa_priors(list(), "cor", "unit_vector", 2:3, 9L)
  expect_identical(priors[[1L]]$h2, .h2_prior_lkj(9L, 2L))
  expect_identical(priors[[2L]]$h2, .h2_prior_lkj(9L, 3L))
  # A prior given by the user applies to every number of factors.
  priors <- .befa_priors(list(h2 = c(1, 1)), "cor", "unit_vector", 2:3, 9L)
  expect_identical(priors[[2L]]$h2, c(1, 1))
})

test_that("the models share the saturated and null models", {
  fits <- lapply(1:3, function(M) .fake_befa("cov", N = 30L, n_iter = 5L, M = M))
  draws <- .null_draws(5L, raw = FALSE)
  n_fits <- 0L
  n_saturated <- 0L
  real <- .saturated_moments
  local_mocked_bindings(
    .run_stan = function(...) {
      n_fits <<- n_fits + 1L
      list(draws = draws)
    },
    .saturated_moments = function(...) {
      n_saturated <<- n_saturated + 1L
      real(...)
    }
  )
  loo_config <- validate_loo_args(NULL)
  shared <- new.env(parent = emptyenv())
  together <- lapply(fits, .fit_measures, loo_config = loo_config, shared = shared)
  expect_identical(n_fits, 1L)
  expect_identical(n_saturated, 1L)

  # The shared models give the indices of a fit on its own.
  alone <- .fit_measures(fits[[3L]], loo_config)
  expect_identical(n_fits, 2L)
  expect_equal(together[[3L]]$posterior_fit, alone$posterior_fit)
  expect_equal(together[[3L]]$information, alone$information)
})

test_that("befa() returns one fit per number of factors, compared", {
  calls <- .mock_fits()
  fits <- .fake_list(c(3, 1, 2), seed = 11)
  expect_s3_class(fits, c("befa_list", "list"), exact = TRUE)
  expect_named(fits, c("nf1", "nf2", "nf3"))
  expect_true(all(vapply(fits, inherits, TRUE, "befa")))
  expect_identical(fits$nf2$n_factors, 2L)
  expect_identical(attr(fits, "n_factors"), 1:3)
  expect_identical(attr(fits, "call")$n_factors, quote(n_factors))
  expect_length(attr(fits, "failed"), 0L)

  # Every model has the seed of the call, the same environment for what the
  # models share, and a call for its own number of factors.
  expect_identical(vapply(calls$run, function(r) r$settings$seed, 0), rep(11, 3))
  expect_identical(calls$shared[[1L]], calls$shared[[3L]])
  expect_identical(lapply(calls$run, function(r) r$call$n_factors), list(1L, 2L, 3L))

  cmp <- attr(fits, "comparison")
  expect_s3_class(cmp, "befa_comparison")
  expect_setequal(cmp$criteria, c("fit", "loo", "ic"))
  expect_null(cmp$evidence)
  expect_setequal(rownames(cmp$fit), names(fits))
  expect_identical(befa_compare(fits), cmp)
})

test_that("compare chooses whether and how the models are compared", {
  .mock_fits()
  expect_null(attr(.fake_list(1:2, compare = FALSE), "comparison"))

  one <- .fake_list(1:2, compare_criteria = "ic")
  expect_identical(attr(one, "comparison")$criteria, "ic")

  local_mocked_bindings(befa_bridge_sampler = function(object, ...) {
    structure(list(logml = -100 - 10 / object$n_factors, mcse_logml = 0.05),
              class = "bridge")
  })
  fits <- .fake_list(1:2, compare_criteria = c("fit", "evidence"))
  expect_setequal(attr(fits, "comparison")$criteria, c("fit", "evidence"))
  overview <- summary(fits)$overview
  expect_true("P(M|y)" %in% names(overview))
  expect_gt(overview$`P(M|y)`[2L], overview$`P(M|y)`[1L])
  expect_match(.shown(summary(fits)), "P(M|y)", fixed = TRUE)
})

test_that("a model that fails does not stop the others", {
  .mock_fits(fail = 2L)
  expect_warning(fits <- .fake_list(1:3), "Models not fitted: nf2 \\(did not converge\\)")
  expect_named(fits, c("nf1", "nf3"))
  expect_identical(attr(fits, "n_factors"), c(1L, 3L))
  expect_identical(attr(fits, "failed"), c(nf2 = "did not converge"))
  expect_setequal(rownames(attr(fits, "comparison")$fit), c("nf1", "nf3"))
  expect_match(.shown(fits), "Not fitted: nf2 (did not converge)", fixed = TRUE)
  expect_match(.shown(summary(fits)), "Not fitted: nf2", fixed = TRUE)
})

test_that("befa() stops when no model can be fitted", {
  .mock_fits(fail = 1:2)
  expect_error(.fake_list(1:2), "No model could be fitted.*nf1: did not converge")
})

test_that("one model left is not compared, and a failed comparison only warns", {
  .mock_fits(fail = 2L)
  expect_warning(fits <- .fake_list(1:2), "Models not fitted")
  expect_length(fits, 1L)
  expect_null(attr(fits, "comparison"))
})

test_that("a comparison that fails keeps the fits", {
  .mock_fits()
  local_mocked_bindings(befa_compare = function(...) stop("no loo"))
  expect_warning(fits <- .fake_list(1:2), "could not be compared: no loo")
  expect_length(fits, 2L)
  expect_null(attr(fits, "comparison"))
})

test_that("warnings name the model that raised them", {
  expect_warning(res <- .with_model_label({
    warning("short chains", call. = FALSE)
    "fit"
  }, "nf2"), "^nf2: short chains$")
  expect_identical(res, "fit")
})

test_that("each model is announced before its progress bar", {
  calls <- .mock_fits()
  # This function returns the messages of a verbose fit as one string.
  messages_of <- function(...) {
    msgs <- character(0)
    withCallingHandlers(
      befa(.simulated_data(N = 60L, seed = 7), 1:2, verbose = TRUE, ...),
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
      })
    paste(msgs, collapse = "")
  }

  msgs <- messages_of()
  expect_identical(vapply(calls$run, `[[`, "", "title"),
                   c("Estimating model 1 of 2 (1 factor)",
                     "Estimating model 2 of 2 (2 factors)"))
  expect_match(msgs, "Comparing the models")
  # The title goes on top of the bar, so it is not repeated.
  expect_no_match(msgs, "Estimating model")

  # Without a progress bar, the title is shown on its own.
  msgs <- messages_of(refresh = 0)
  expect_match(msgs, "Estimating model 1 of 2 (1 factor)", fixed = TRUE)
  expect_match(msgs, "Estimating model 2 of 2 (2 factors)", fixed = TRUE)
})

test_that("the header of the fit shows its title", {
  settings <- .sampler_settings("mcmc", list(chains = 2, seed = 1))
  expect_message(.stan_header("cmdstanr", settings, "Estimating model 1 of 3"),
                 "Estimating model 1 of 3", fixed = TRUE)
  expect_message(.stan_header("cmdstanr", settings),
                 "Bayesian exploratory factor analysis", fixed = TRUE)
})

test_that("the null model is announced only while it is pending", {
  object <- .fake_befa()
  local_mocked_bindings(.fit_measures = function(...) NULL)
  shared <- new.env(parent = emptyenv())
  expect_message(.befa_downstream(object, FALSE, FALSE, TRUE, 1L, TRUE, shared),
                 "Computing the fit indices and the null model.", fixed = TRUE)
  shared$null <- list(df = 1)
  expect_message(.befa_downstream(object, FALSE, FALSE, TRUE, 1L, TRUE, shared),
                 "Computing the fit indices.", fixed = TRUE)
})

test_that("summary statistics warn once and compare only the evidence", {
  n_warn <- 0L
  prepared <- list(summary = TRUE)
  local_mocked_bindings(.befa_fit_one = function(inputs, run, shared = NULL) {
    expect_false(run$factor_scores)
    expect_false(run$compute_fit_indices)
    .fake_befa(M = inputs$n_factors)
  })
  run <- list(prepared = prepared, factor_scores = TRUE, compute_fit_indices = TRUE,
              verbose = FALSE, call = quote(befa()), settings = list(seed = 1))
  inputs <- list(list(n_factors = 1L), list(n_factors = 2L))
  withCallingHandlers(
    fits <- .befa_fit_list(inputs, run, c("fit", "loo", "ic")),
    warning = function(w) {
      n_warn <<- n_warn + 1L
      expect_match(conditionMessage(w), "need the raw `data`")
      invokeRestart("muffleWarning")
    })
  expect_identical(n_warn, 2L)
  expect_null(attr(fits, "comparison"))

  # Without fit indices or ELPD, no model is chosen to be shown.
  s <- summary(fits)
  expect_length(s$shown, 0L)
  expect_named(s$overview, c("Model", "Rhat", "ESS"))
  expect_match(.shown(s), "No model has an ELPD to choose from", fixed = TRUE)
  expect_no_match(.shown(s), "Standardized Factor Loadings", fixed = TRUE)
})

test_that("without fit indices the overview reads the comparison", {
  local_mocked_bindings(.befa_fit_one = function(inputs, run, shared = NULL) {
    .fake_befa(M = inputs$n_factors)
  })
  run <- list(prepared = list(summary = FALSE), factor_scores = FALSE,
              compute_fit_indices = FALSE, verbose = FALSE, call = quote(befa()),
              settings = list(seed = 1))
  inputs <- list(list(n_factors = 1L), list(n_factors = 2L))
  fits <- .befa_fit_list(inputs, run, c("fit", "loo", "ic"))
  cmp <- attr(fits, "comparison")
  expect_setequal(cmp$criteria, c("loo", "ic"))
  overview <- summary(fits)$overview
  expect_false("BRMSEA" %in% names(overview))
  expect_equal(overview$ELPD, unname(cmp$loo$table[c("nf1", "nf2"), "ELPD"]))
})

test_that("the overview holds one row per model and fits in 80 columns", {
  .mock_fits()
  fits <- .fake_list()
  s <- summary(fits)
  expect_s3_class(s, "summary.befa_list")
  o <- s$overview
  expect_named(o, c("Model", "Rhat", "ESS", "BRMSEA", "BCFI", "ELPD", "WAIC"))
  expect_identical(o$Model, c("nf1", "nf2", "nf3"))
  fm <- fits$nf2$fit_indices
  expect_equal(o$BRMSEA[2L], mean(fm$posterior_fit$BRMSEA))
  expect_equal(o$ELPD[2L], fm$information["ELPD", "Estimate"])
  expect_equal(o$WAIC[2L], fm$information["WAIC", "Estimate"])

  # The widest table has every column and the numbers of a large data set.
  wide <- data.frame(Model = "nf10", Rhat = 1.003, ESS = 12345, Div = 123,
                     BRMSEA = 0.123, BCFI = 0.987, ELPD = -123456.78,
                     WAIC = 246913.56, `P(M|y)` = 1, check.names = FALSE)
  lines <- utils::capture.output(width <- .print_overview(wide, 2, "Table"))
  expect_lte(width, 80L)
  expect_true(all(nchar(trimws(lines, "right")) <= 80L))
})

test_that("approximations report their Pareto k", {
  .mock_fits()
  fits <- .fake_list(1:2)
  for (k in names(fits)) {
    fits[[k]]$diagnostics$max_rhat <- NA_real_
    fits[[k]]$diagnostics$min_ess_bulk <- NA_real_
    fits[[k]]$diagnostics$min_ess_tail <- NA_real_
    fits[[k]]$diagnostics$pareto_k <- 0.4
  }
  o <- summary(fits)$overview
  expect_true("Pareto k" %in% names(o))
  expect_false("Rhat" %in% names(o))
  expect_match(.shown(fits), "Pareto k of the approximation", fixed = TRUE)
})

test_that("the summary shows only the model with the highest ELPD by default", {
  .mock_fits()
  fits <- .fake_list()
  s <- summary(fits)
  best <- s$overview$Model[which.max(s$overview$ELPD)]
  expect_identical(s$shown, best)
  expect_named(s$models, best)
  expect_s3_class(s$models[[best]], "summary.befa")

  out <- .shown(s)
  for (text in c("Bayesian EFA with 1, 2, 3 factors",
                 "Data model: cor; prior: unit_vector; N = 60",
                 "Table 1. Overview of the Models",
                 "Table 2. Bayesian Fit Measures",
                 "Table 3. LOO Comparison",
                 "Table 4. Information Criteria",
                 sprintf("Table 5. Standardized Factor Loadings (%s)", best),
                 sprintf("Table 6. Variance Explained (%s)", best),
                 sprintf("Only %s is shown, the model with the highest ELPD.", best),
                 "the posterior may be unreliable in")) {
    expect_match(out, text, fixed = TRUE)
  }
  expect_no_match(out, "Table 7", fixed = TRUE)
  expect_output(res <- print(s), "Table 1")
  expect_identical(res, s)
})

test_that("show chooses the models and the options reach their tables", {
  .mock_fits()
  fits <- .fake_list()
  expect_identical(summary(fits, show = "all")$shown, c("nf1", "nf2", "nf3"))
  expect_identical(summary(fits, show = c(3, 1))$shown, c("nf1", "nf3"))
  expect_identical(summary(fits, show = "nf2")$shown, "nf2")
  expect_error(summary(fits, show = 7), "fitted \\(1, 2, 3\\)")
  expect_error(summary(fits, show = list(1)), "`show`")
  expect_error(summary(fits, show = character(0)), "`show`")

  s <- summary(fits, show = 2, cutoff = 0.3, signif_stars = TRUE)
  expect_null(s$best)
  out <- .shown(s)
  expect_match(out, "Standardized Factor Loadings (nf2)", fixed = TRUE)
  expect_match(out, "are hidden", fixed = TRUE)
  expect_no_match(out, "Only", fixed = TRUE)
  expect_no_match(.shown(s, cutoff = 0, signif_stars = FALSE), "hidden", fixed = TRUE)

  all <- .shown(summary(fits, show = "all"))
  expect_match(all, "Table 10. Variance Explained (nf3)", fixed = TRUE)
})

test_that("models with fewer factors close to the best one are pointed out", {
  tab <- cbind(elpd_diff = c(0, -3, -40), se_diff = c(0, 2, 5))
  rownames(tab) <- c("nf3", "nf2", "nf1")
  cmp <- list(loo = list(table = tab))
  expect_identical(.simpler_models(cmp, "nf3"), "nf2")
  expect_identical(.simpler_models(cmp, "nf2"), character(0))
  expect_identical(.simpler_models(list(), "nf3"), character(0))
  expect_identical(.simpler_models(cmp, character(0)), character(0))

  tab[, "se_diff"] <- c(0, 2, 30)
  cmp$loo$table <- tab
  expect_identical(.simpler_models(cmp, "nf3"), c("nf1", "nf2"))

  best <- list(model = "nf3", criterion = "ELPD")
  note <- paste(.shown_note(best, "nf2"), collapse = " ")
  expect_match(note, "nf2 has fewer factors and an ELPD within 2 SE of it.", fixed = TRUE)
  note <- paste(.shown_note(best, c("nf1", "nf2")), collapse = " ")
  expect_match(note, "nf1, nf2 have fewer factors", fixed = TRUE)
  note <- paste(.shown_note(list(model = "nf2", criterion = "P(M|y)"), character(0)),
                collapse = " ")
  expect_match(note, "highest posterior probability", fixed = TRUE)
  expect_no_match(note, "fewer", fixed = TRUE)
})

test_that("the best model falls back on the posterior probability", {
  o <- data.frame(Model = c("nf1", "nf2"), `P(M|y)` = c(0.1, 0.9), check.names = FALSE)
  expect_identical(.best_model(o), list(model = "nf2", criterion = "P(M|y)"))
  o$ELPD <- c(-10, -20)
  expect_identical(.best_model(o), list(model = "nf1", criterion = "ELPD"))
})

test_that("print() shows the overview only", {
  .mock_fits()
  fits <- .fake_list()
  expect_output(res <- print(fits), "Overview of the Models")
  expect_identical(res, fits)
  out <- .shown(fits)
  expect_match(out, "x$nf1 is the fit of that model.", fixed = TRUE)
  expect_no_match(out, "Standardized Factor Loadings", fixed = TRUE)
})

test_that("the functions of a fit work on the list or point to its models", {
  .mock_fits()
  fits <- .fake_list(1:2)
  fm <- befa_fit_measures(fits)
  expect_named(fm, c("nf1", "nf2"))
  expect_s3_class(fm$nf2, "befa_fitmeasures")
  rel <- befa_reliability(fits)
  expect_named(rel, c("nf1", "nf2"))
  expect_s3_class(rel$nf1, "befa_reliability")
  expect_error(extract_posterior_draws(fits), "x$nf1", fixed = TRUE)
  expect_error(befa_factor_scores(fits), "several fits")

  # The list is compared again when other arguments are given, or a part of it.
  ic <- befa_compare(fits, criteria = "ic")
  expect_identical(ic$criteria, "ic")
  sub <- befa_compare(fits[1:2], criteria = "fit")
  expect_identical(sub$models, c("nf1", "nf2"))
})

test_that("befa() fits one model per number of factors with stanli", {
  skip_on_cran()
  skip_if_not_installed("stanli", "0.25.0")
  skip_if_not(stanli::stanli_available(), "the stanli runtime is not installed")
  withr::local_options(BayesianEFA.cache_dir = withr::local_tempdir())
  Y <- .simulated_data(N = 300)
  args <- list(backend = "stanli", chains = 2, iter_warmup = 200,
               iter_sampling = 150, seed = 5, verbose = FALSE)
  fits <- suppressWarnings(do.call(befa, c(list(Y, 1:2, model = "cov"), args)))
  alone <- suppressWarnings(do.call(befa, c(list(Y, 2, model = "cov"), args)))
  expect_s3_class(fits, "befa_list")
  expect_equal(fits$nf2$draws, alone$draws)
  # The shared null model gives the indices of the model fitted on its own.
  expect_equal(fits$nf2$fit_indices$posterior_fit, alone$fit_indices$posterior_fit)
  expect_setequal(attr(fits, "comparison")$criteria, c("fit", "loo", "ic"))
  expect_gt(summary(fits)$overview$ELPD[2L], summary(fits)$overview$ELPD[1L])
})
