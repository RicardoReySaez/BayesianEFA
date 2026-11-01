# Comparison of fitted models

# This function returns a fake fit to the same data as .fake_befa() whose
# standardized loadings are shrunk by `f`, so it fits worse.
.shrunk_befa <- function(f = 0.6) {
  object <- .fake_befa()
  L <- .stacked_draws(object$draws, "Lambda_std") * f
  h2 <- .item_sums(L^2, 5L, 2L)
  dims <- dim(object$draws)
  object$draws <- .as_chains(.named_draws(list(Lambda_std = L, Psi_std = 1 - h2, h2 = h2),
                                          5L, 2L), dims[1L], dims[2L])
  object
}

# This function returns fake bridge results with the log evidences and Monte
# Carlo errors of their repetitions.
.fake_bridge <- function(logml, mcse = rep(0.05, length(logml))) {
  structure(list(logml = logml, mcse_logml = mcse),
            class = if (length(logml) > 1L) "bridge_list" else "bridge")
}

test_that("the LOO table follows loo and is ordered by ELPD", {
  good <- .fake_befa()
  bad <- .shrunk_befa()
  cmp <- befa_compare(bad = bad, good = good, criteria = c("loo", "ic"), seed = 1)
  expect_s3_class(cmp, "befa_comparison")
  expect_identical(cmp$models, c("bad", "good"))

  tab <- cmp$loo$table
  expect_identical(rownames(tab), c("good", "bad"))
  expect_identical(unname(tab["good", "elpd_diff"]), 0)
  expect_lt(tab["bad", "elpd_diff"], 0)
  ref <- as.data.frame(loo::loo_compare(cmp$loo$loos))
  ref_bad <- if ("model" %in% names(ref)) ref[ref$model == "bad", ] else ref["bad", ]
  expect_equal(unname(tab["bad", c("elpd_diff", "se_diff")]),
               c(ref_bad$elpd_diff, ref_bad$se_diff))
  expect_equal(unname(tab["bad", "p_worse"]),
               stats::pnorm(0, ref_bad$elpd_diff, ref_bad$se_diff))
  expect_true(is.na(tab["good", "p_worse"]))
  expect_equal(sum(tab[, "w_stacking"]), 1, tolerance = 1e-6)
  expect_equal(sum(tab[, "w_pseudobma"]), 1, tolerance = 1e-6)
  expect_gt(tab["good", "w_pseudobma"], tab["bad", "w_pseudobma"])

  # Stored fit indices are reused, and the same seed gives the same weights.
  good$fit_indices <- befa_fit_measures(good)
  again <- befa_compare(bad = bad, good = good, criteria = c("loo", "ic"), seed = 1)
  expect_equal(again$loo$table, cmp$loo$table)
})

test_that("the criteria come with Akaike weights, ordered by WAIC", {
  good <- .fake_befa()
  bad <- .shrunk_befa()
  cmp <- befa_compare(bad = bad, good = good, criteria = "ic")
  tab <- cmp$ic
  expect_identical(rownames(tab), c("good", "bad"))
  expect_identical(colnames(tab), c("WAIC", "w_WAIC", "DIC", "w_DIC", "DIC_i",
                                    "w_DIC_i", "BPIC", "w_BPIC"))
  fm <- befa_fit_measures(good)
  expect_equal(tab["good", "WAIC"], fm$information["WAIC", "Estimate"])
  waic <- tab[, "WAIC"]
  expect_equal(unname(tab[, "w_WAIC"]),
               unname(exp(-(waic - min(waic)) / 2) / sum(exp(-(waic - min(waic)) / 2))))
  expect_equal(.akaike_weights(c(10, 12)), c(1, exp(-1)) / (1 + exp(-1)))
  expect_true(all(tab["good", c("WAIC", "DIC", "DIC_i", "BPIC")] <
                    tab["bad", c("WAIC", "DIC", "DIC_i", "BPIC")]))
})

test_that("the fit indices are compared by their posterior means", {
  # The posterior predictive p-value draws replicated data, so the stored
  # fit measures are the reference.
  good <- .fake_befa()
  fm <- befa_fit_measures(good)
  good$fit_indices <- fm
  bad <- .shrunk_befa()
  cmp <- befa_compare(bad = bad, good = good, criteria = "fit")
  tab <- cmp$fit
  expect_identical(rownames(tab), c("good", "bad"))
  expect_identical(colnames(tab), c("Chi2_ppp", "BRMSEA", "BGamma", "Adj_BGamma",
                                    "BMc", "SRMR", "BCFI", "BTLI"))
  expect_equal(tab["good", "BRMSEA"], mean(.stacked_draws(fm$fit_draws, "BRMSEA")))
  expect_equal(tab["good", "Chi2_ppp"], fm$details$chi2_ppp)
  expect_lte(tab["good", "BRMSEA"], tab["bad", "BRMSEA"])
})

test_that("the evidence gives posterior model probabilities", {
  local_mocked_bindings(befa_bridge_sampler = function(object, ...) {
    if (isTRUE(object$fail)) stop("the prior is not invariant")
    if (isTRUE(object$weak)) .fake_bridge(-110, 0.08) else .fake_bridge(-100, 0.04)
  })
  bad <- .fake_befa()
  bad$weak <- TRUE
  good <- .fake_befa()
  cmp <- befa_compare(bad, good, criteria = "evidence", seed = 2)
  tab <- cmp$evidence$table
  expect_identical(rownames(tab), c("Model 2", "Model 1"))
  expect_equal(unname(tab[, "log_ML"]), c(-100, -110))
  expect_equal(unname(tab[, "MC_SE"]), c(0.04, 0.08))
  expect_equal(unname(tab[, "Estimate"]), c(1, exp(-10)) / (1 + exp(-10)))
  expect_true(all(tab[, "CI_Low"] <= tab[, "Estimate"] & tab[, "Estimate"] <= tab[, "CI_High"]))

  # Prior model probabilities shift the posterior ones.
  tilted <- befa_compare(bad, good, criteria = "evidence", model_prior = c(1e6, 1))
  expect_gt(tilted$evidence$table["Model 1", "Estimate"], tab["Model 1", "Estimate"])
  expect_equal(unname(tilted$evidence$model_prior), c(1e6, 1) / (1 + 1e6))

  # With repetitions, the medians of the estimates and of their errors
  local_mocked_bindings(befa_bridge_sampler = function(object, ...) {
    .fake_bridge(c(-100, -100.2, -99.9), c(0.04, 0.06, 0.05))
  })
  reps <- befa_compare(bad, good, criteria = "evidence", repetitions = 3)
  expect_equal(unname(reps$evidence$table[, "log_ML"]), c(-100, -100))
  expect_equal(unname(reps$evidence$table[, "MC_SE"]), c(0.05, 0.05))

  # An unreliable Monte Carlo error is noted.
  local_mocked_bindings(befa_bridge_sampler = function(object, ...) .fake_bridge(-100, 0.6))
  noisy <- befa_compare(bad, good, criteria = "evidence")
  expect_match(noisy$notes, "above 0.3")

  # A model without evidence is noted and leaves no probabilities.
  local_mocked_bindings(befa_bridge_sampler = function(object, ...) {
    if (isTRUE(object$fail)) stop("the prior is not invariant")
    .fake_bridge(-100)
  })
  bad$fail <- TRUE
  failed <- befa_compare(bad, good, criteria = "evidence")
  expect_null(failed$evidence$table)
  expect_match(failed$notes, "No evidence for Model 1: the prior is not invariant")
  expect_output(print(failed), "No evidence for Model 1", fixed = TRUE)

  # Without Monte Carlo errors, the probabilities have no interval.
  local_mocked_bindings(befa_bridge_sampler = function(object, ...) list(logml = -100))
  plain <- befa_compare(good, good, criteria = "evidence")
  expect_true(all(is.na(plain$evidence$table[, "CI_Low"])))
})

test_that("the comparison prints four tables within their width", {
  local_mocked_bindings(befa_bridge_sampler = function(object, ...) {
    .fake_bridge(-100 - 10 * isTRUE(object$weak), 0.05)
  })
  good <- .fake_befa()
  bad <- .shrunk_befa()
  bad$weak <- TRUE
  cmp <- befa_compare(bad = bad, a_very_long_model_name = good, seed = 3)
  out <- utils::capture.output(res <- print(cmp))
  expect_identical(res, cmp)
  text <- paste(out, collapse = "\n")
  for (s in c("Table 1. Bayesian Fit Measures", "Table 2. LOO Comparison",
              "Table 3. Information Criteria", "Table 4. Model Evidence",
              "pseudo-BMA+", "Akaike weights", "Micaletto and Vehtari",
              "a_ver…name")) {
    expect_match(text, s, fixed = TRUE)
  }
  # Every table lists the best model first.
  first_rows <- out[which(grepl("^—", out)) + 1L]
  expect_true(all(startsWith(first_rows, "a_ver…name")))
  # Every line of a table and of its note, which ends at the next blank
  # line, fits within the table's width.
  rules <- which(grepl("^‗", out))
  for (k in seq(1, length(rules), by = 2)) {
    width <- nchar(trimws(out[rules[k]], "right"))
    end <- rules[k + 1L]
    while (end < length(out) && nzchar(trimws(out[end + 1L]))) end <- end + 1L
    expect_true(all(nchar(trimws(out[rules[k]:end], "right")) <= width))
  }
})

test_that("only fits to the same data and data model are compared", {
  good <- .fake_befa()
  expect_error(befa_compare(good), "two or more fits")
  expect_error(befa_compare(good, list()), "two or more fits")
  expect_error(befa_compare(good, .fake_befa("cov")), "share the data model")
  expect_error(befa_compare(good, .fake_befa(N = 50L)), "same data")
  expect_error(befa_compare(good, good, model_prior = c(1, -1)), "positive numbers")
  expect_identical(befa_compare(list(good, x = good), criteria = "ic")$models,
                   c("Model 1", "x"))

  summ <- good
  summ$data$summary <- TRUE
  summ$data$S <- diag(5)
  local_mocked_bindings(befa_bridge_sampler = function(object, ...) .fake_bridge(-1))
  only <- befa_compare(summ, summ)
  expect_identical(only$criteria, "evidence")
  expect_match(only$notes[1], "summary statistics")
})

test_that("befa_compare() prefers the true number of factors", {
  skip_on_cran()
  skip_if_not_installed("bridgesampling", "1.2.1")
  skip_if_not_installed("stanli", "0.25.0")
  skip_if_not(stanli::stanli_available(), "the stanli runtime is not installed")
  withr::local_options(BayesianEFA.cache_dir = withr::local_tempdir())
  Y <- .simulated_data(N = 300)
  fit <- function(M) suppressWarnings(befa(Y, M, backend = "stanli", chains = 2,
                                           iter_warmup = 500, iter_sampling = 1500,
                                           seed = 6, verbose = FALSE,
                                           compute_reliability = FALSE))
  cmp <- befa_compare(one = fit(1L), two = fit(2L), seed = 7)
  expect_identical(rownames(cmp$fit)[1], "two")
  expect_identical(rownames(cmp$loo$table)[1], "two")
  expect_gt(cmp$loo$table["two", "w_stacking"], 0.5)
  expect_identical(rownames(cmp$ic)[1], "two")
  expect_gt(cmp$ic["two", "w_WAIC"], 0.9)
  expect_identical(rownames(cmp$evidence$table)[1], "two")
  expect_gt(cmp$evidence$table["two", "Estimate"], 0.9)
  expect_true(all(cmp$evidence$table[, "MC_SE"] < 0.3))
})
