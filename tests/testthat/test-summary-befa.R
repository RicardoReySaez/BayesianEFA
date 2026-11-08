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
  draws <- .stacked_draws(object$draws, "Lambda_std")
  expect_equal(s$variance$per_factor, colSums(matrix(colMeans(draws^2), 5, 2)) / 5)
  expect_equal(s$variance$total_explained, sum(s$variance$per_factor))
  # The total of the variance table is the sum of the communalities.
  expect_equal(s$variance$total_explained, sum(s$estimates$h2) / 5)
  variance <- s$tables$table2_variance
  expect_identical(colnames(variance), c("F1", "F2", "Total"))
  expect_equal(variance["Proportion of total", "F1"] +
                 variance["Proportion of total", "F2"], 1)
  expect_equal(variance["Cumulative variance", "F2"], s$variance$total_explained)
  expect_true(is.na(variance["Cumulative variance", "Total"]))
  expect_true(all(s$sig_matrix[cbind(1:5, c(1, 1, 1, 2, 2))]))
  expect_identical(s$header$N, 60L)
  expect_identical(names(s$tables),
                   c("table1_loadings", "table2_variance", "table3_fit_measures",
                     "table4_information", "table5_reliability"))
  expect_identical(rownames(s$tables$table4_information)[c(1, 8)], c("ELPD", "DIC_i"))
  expect_identical(rownames(s$tables$table5_reliability), c("F1", "F2", "Total"))
  total <- unlist(s$tables$table5_reliability["Total", ])
  expect_equal(unname(total[1:2]), c(object$reliability$total$Estimate,
                                     object$reliability$total$SD))
  expect_identical(s$estimates$Variable, paste0("Item_", 1:5))
  expect_true(s$print_options$sort)
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
                 "Table 2. Variance Explained",
                 "Sum of squared loadings",
                 "Table 3. Bayesian Fit Measures",
                 "Table 4. Information Criteria",
                 "SEs across persons (jackknife for DIC_i and p_V)",
                 "Table 5. Factor Reliability (Omega)",
                 "Total is the omega of the full scale",
                 "varimax rotation applied; MCMC with stanli",
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
  expect_identical(names(s$tables), c("table1_loadings", "table2_variance"))
  out <- .printed(s)
  expect_match(out, "Table 2.", fixed = TRUE)
  expect_no_match(out, "Table 3.", fixed = TRUE)
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

test_that("long item names are shortened in the middle and stay unique", {
  expect_identical(.short_names(c("x1", "Variable")), c("x1", "Variable"))
  expect_identical(.short_names("neuroticism_12"), "neur\u2026_12")
  # Two names that would become equal end with their item number instead.
  short <- .short_names(c("speeded_addition", "cubes", "speeded_discrimination"))
  expect_identical(short, c("speede\u20261", "cubes", "speede\u20263"))
  expect_true(all(nchar(short) <= 8L))
  expect_identical(.short_names("neuroticism_12", 6L), "neu\u202612")
})

test_that("the loadings table keeps its width whatever the item names", {
  object <- .fake_befa()
  # This function returns the width of the first table that print() draws.
  .table_width <- function(s, ...) {
    out <- utils::capture.output(print(s, ...))
    nchar(trimws(out[grepl("^\u2017", out)][1L], "right"))
  }
  short <- summary(object)
  object$data$item_names <- c("visual_perception", "cubes", "lozenges",
                              "paragraph_comprehension", "sentence_completion")
  long <- summary(object)
  expect_identical(long$estimates$Variable, object$data$item_names)
  expect_identical(.table_width(long), .table_width(short))
  out <- paste(utils::capture.output(print(long)), collapse = "\n")
  expect_match(out, "visu\u2026ion", fixed = TRUE)
  expect_match(out, "lozenges", fixed = TRUE)
  wide <- .table_width(long, name_width = 12)
  expect_identical(wide, .table_width(short) + 4L)
  # Every note line fits within the table.
  lines <- utils::capture.output(print(long))
  expect_true(all(nchar(trimws(lines, "right")) <= .table_width(long)))
  expect_error(summary(object, name_width = 2), "`name_width`")
})

test_that("the items are sorted by their main factor unless asked not to", {
  object <- .fake_befa()
  # The fake loadings put items 1-3 on F1 and 4-5 on F2, where item 5 loads
  # more than item 4, so it comes first.
  object$data$item_names <- paste0("q", 1:5)
  rows <- function(s, ...) {
    out <- utils::capture.output(print(s, ...))
    sub(" .*", "", out[grepl("^q[0-9]", out)])
  }
  s <- summary(object)
  expect_identical(rows(s), c("q1", "q2", "q3", "q5", "q4"))
  expect_identical(rows(s, sort = FALSE), paste0("q", 1:5))
})

test_that("summaries of fits from earlier versions still print", {
  object <- .fake_befa()
  object$data$item_names <- NULL
  s <- summary(object)
  expect_identical(s$estimates$Variable, paste0("Item_", 1:5))
  s$print_options$name_width <- NULL
  expect_output(print(s), "Item_1", fixed = TRUE)
})
