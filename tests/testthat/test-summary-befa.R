# Summary of a fit and its printed tables

# This function returns a fake fit with its omega coefficients and fit
# indices, as befa() stores them.
.fake_full_befa <- function(...) {
  object <- .fake_befa(...)
  object$reliability <- befa_reliability(object)
  object$fit_indices <- befa_fit_measures(object)
  object
}

# This function returns what print() shows as one line, so wrapped notes
# can be matched.
.printed <- function(x, ...) {
  gsub("\\s+", " ", paste(utils::capture.output(print(x, ...)), collapse = " "))
}

test_that("the summary holds the standardized solution of the fit", {
  object <- .fake_full_befa()
  s <- summary(object)
  expect_s3_class(s, "summary.befa")
  expect_named(s$estimates, c("Variable", "F1", "F2", "h2", "u2", "Rhat",
                              "EssBulk", "EssTail"))
  L <- matrix(colMeans(.stacked_draws(object$draws, "Lambda_std")), 5, 2)
  expect_equal(as.matrix(s$estimates[, c("F1", "F2")]), L, ignore_attr = TRUE)
  expect_equal(s$estimates$h2 + s$estimates$u2, rep(1, 5))
  expect_equal(s$variance$per_factor, colSums(L^2) / 5)
  expect_equal(s$variance$total_explained, sum(s$variance$per_factor))
  expect_true(all(s$sig_matrix[cbind(1:5, c(1, 1, 1, 2, 2))]))
  expect_identical(s$header$N, 60L)
  expect_identical(names(s$tables), c("table1_loadings", "table2_fit_measures",
                                      "table3_reliability"))
  expect_identical(rownames(s$tables$table3_reliability), c("F1", "F2"))
})

test_that("the summary uses the standardized metric of every model", {
  cor <- summary(.fake_befa("cor"))
  cov <- summary(.fake_befa("cov", lambda_prior = "normal"))
  expect_equal(cov$estimates[, 2:5], cor$estimates[, 2:5])
})

test_that("the summary prints its tables with their options", {
  object <- .fake_full_befa()
  out <- .printed(summary(object))
  for (text in c("Table 1. Standardized Factor Loadings",
                 "Table 2. Bayesian Fit Measures",
                 "Table 3. Factor Reliability (Coefficient Omega)",
                 "Full Scale Omega Total", "varimax rotation applied; MCMC with stanli",
                 "the posterior may be unreliable")) {
    expect_match(out, text, fixed = TRUE)
  }

  s <- summary(object, cutoff = 0.3, sort = TRUE, signif_stars = TRUE)
  expect_output(res <- print(s), "Table 1")
  expect_identical(res, s)
  out <- .printed(s, digits = 3)
  expect_match(out, "Loadings with absolute values < 0.30 are hidden.", fixed = TRUE)
  expect_match(out, "(*) 95% Credible Interval excludes 0.", fixed = TRUE)
  expect_match(out, "0.\\d{3}\\*")
  out <- .printed(s, cutoff = 0, signif_stars = FALSE)
  expect_no_match(out, "hidden", fixed = TRUE)
})

test_that("a fit without fit indices or reliability prints its loadings only", {
  object <- .fake_befa()
  object$diagnostics$issues <- character(0)
  s <- summary(object)
  expect_identical(names(s$tables), "table1_loadings")
  out <- .printed(s)
  expect_match(out, "Table 1.", fixed = TRUE)
  expect_no_match(out, "Table 2.", fixed = TRUE)
  expect_no_match(out, "unreliable", fixed = TRUE)
})

test_that("the table helpers sort, format and split their values", {
  expect_identical(.format_numeric(c(1.234, NA, 2)), c("1.23", "", "2.00"))
  expect_identical(.format_numeric(list("a", ""), 1), c("a", ""))

  df <- data.frame(Variable = c("a", "b", "c"), F1 = c(0.1, 0.9, 0.2),
                   F2 = c(0.8, 0.1, 0.05), EssBulk = c(400.4, 500, 600))
  sig <- matrix(c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE), 3)
  res <- .format_loading_table(df, sig, c("F1", "F2"), cutoff = 0.15,
                               sort = TRUE, signif_stars = TRUE, digits = 2)
  expect_identical(res$df$Variable, c("b", "c", "a"))
  expect_identical(res$df$F1, c("0.90*", "0.20*", ""))
  expect_identical(res$df$F2, c("", "", "0.80*"))
  expect_identical(res$df$EssBulk, c("500", "600", "400"))

  rel <- build_reliability_table(befa_reliability(.fake_befa()))
  expect_named(rel$total, c("Estimate", "SD", "CI_Low", "CI_High"))
  expect_identical(rownames(rel$subscales), c("F1", "F2"))
})
