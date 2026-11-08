# summary-befa.R
# Summary and printed tables of a fitted model
# Author: Ricardo Rey-Sáez
# Last modified: 08-11-2026

# Public functions

#' Summarize BEFA Results
#'
#' Summarises the standardized solution of a fitted model: the posterior
#' means of the aligned loadings, their communalities and uniquenesses, the
#' convergence of each item and the variance each factor explains, followed by
#' the fit indices and the omega coefficients when the fit holds them.
#'
#' @param object A `befa` object returned by [befa()].
#' @param probs Numeric vector of length 2 with the quantiles of the credible
#'   intervals (default `c(0.025, 0.975)`).
#' @param cutoff Loadings whose absolute posterior mean is below `cutoff` are
#'   hidden when printing.
#' @param sort Logical. If `TRUE` (default), the items are grouped by the
#'   factor of their largest loading, from the largest loading to the
#'   smallest within each factor, as `psych::fa.sort()` does. The factors
#'   themselves are always ordered by the variance they explain.
#' @param signif_stars Logical. Mark the loadings whose credible interval
#'   excludes zero with an asterisk.
#' @param name_width Width of the column of item names when printing (default
#'   8, so the table keeps its width). Longer names are shortened in the
#'   middle, which keeps their start and end.
#' @param ... Ignored.
#'
#' @return A `summary.befa` object with:
#' * `$estimates`: One row per item with its full name (from the columns of
#'   `data` or of the summary statistics; `Item_1`, ... when they have none),
#'   the standardized loadings, `h2`, `u2`, the largest R-hat and the
#'   smallest bulk and tail ESS of its loadings.
#' * `$sig_matrix`: Which loadings have a credible interval that excludes 0.
#' * `$variance`: The share of the variance of the standardized items that
#'   each factor, and all of them, account for (posterior means).
#' * `$diagnostics`, `$fit_indices`, `$reliability`: Those of the fit.
#' * `$header`, `$probs`, `$print_options` and `$tables`: What
#'   [print.summary.befa()] shows.
#'
#' @examples
#' \dontrun{
#' fit <- befa(HS_data, n_factors = 3, seed = 17)
#' summary(fit, cutoff = 0.3, signif_stars = TRUE)
#' }
#'
#' @method summary befa
#' @export
summary.befa <- function(object, probs = c(0.025, 0.975), cutoff = 0,
                         sort = TRUE, signif_stars = FALSE, name_width = 8L,
                         ...) {
  .check_befa_object(object)
  .check_count(name_width, "name_width", lower = 4)
  J <- object$stan_data$J
  M <- object$n_factors
  factor_names <- paste0("F", seq_len(M))

  draws <- .metric_draws(object, "std")
  info <- .extract_loading_summary(
    posterior::subset_draws(draws, variable = "Lambda_std"), J, M, probs)
  h2 <- colMeans(.stacked_draws(draws, "h2"))

  item_names <- object$data$item_names
  if (is.null(item_names)) item_names <- .item_labels(NULL, J)
  estimates <- data.frame(
    Variable = item_names, as.data.frame(info$lambda),
    h2 = unname(h2), u2 = 1 - unname(h2),
    Rhat = apply(info$rhat, 1, max, na.rm = TRUE),
    EssBulk = apply(info$ess_bulk, 1, min, na.rm = TRUE),
    EssTail = apply(info$ess_tail, 1, min, na.rm = TRUE)
  )
  colnames(estimates) <- c("Variable", factor_names, "h2", "u2", "Rhat",
                           "EssBulk", "EssTail")
  variance <- .variance_table(.stacked_draws(draws, "Lambda_std"), J, M)

  res <- list(
    estimates = estimates, sig_matrix = info$sig_matrix, probs = probs,
    variance = list(
      per_factor = unname(variance["Proportion of variance", seq_len(M)]),
      total_explained = unname(variance["Proportion of variance", "Total"])
    ),
    diagnostics = object$diagnostics, fit_indices = object$fit_indices,
    reliability = object$reliability,
    header = list(model_type = object$model_type, n_factors = M,
                  N = object$data$N, rotation = object$rotation,
                  method = object$method, backend = object$backend),
    print_options = list(cutoff = cutoff, sort = sort,
                         signif_stars = signif_stars,
                         name_width = as.integer(name_width)),
    tables = list(table1_loadings = estimates,
                  table2_variance = as.data.frame(variance))
  )
  if (!is.null(object$fit_indices)) {
    res$tables$table3_fit_measures <- as.data.frame(
      build_fit_summary_table(object$fit_indices, probs))
    res$tables$table4_information <- as.data.frame(
      build_information_table(object$fit_indices, probs))
  }
  if (!is.null(object$reliability)) {
    rel <- build_reliability_table(object$reliability, probs)
    tab <- rbind(rel$subscales, rel$total)
    rownames(tab) <- c(rownames(rel$subscales), "Total")
    res$tables$table5_reliability <- tab
  }
  structure(res, class = "summary.befa")
}

#' Print BEFA Summary
#'
#' Prints the tables of [summary.befa()]: the standardized loadings, the
#' variance the factors explain, the Bayesian fit measures, the information
#' criteria and the omega coefficients.
#'
#' @param x A `summary.befa` object.
#' @param digits Number of decimal places.
#' @param cutoff,sort,signif_stars,name_width Override the options given to
#'   [summary.befa()].
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#' @method print summary.befa
#' @export
print.summary.befa <- function(x, digits = 2, cutoff = NULL, sort = NULL,
                               signif_stars = NULL, name_width = NULL, ...) {
  options <- .loading_options(x, cutoff, sort, signif_stars, name_width)
  loadings <- .print_loadings(x, digits, options, first = 1L)
  n_table <- loadings$next_table

  # Fit measures and information criteria
  if (!is.null(x$tables$table3_fit_measures)) {
    n_table <- .print_fit_tables(x$tables$table3_fit_measures,
                                 x$tables$table4_information, x$probs,
                                 digits, first = n_table)
  }

  # Reliability
  if (!is.null(x$tables$table5_reliability)) {
    tab <- x$tables$table5_reliability
    out <- as.data.frame(lapply(tab, .format_numeric, digits = digits))
    out <- cbind(Factor = rownames(tab), out)
    cat("\n")
    width <- .print_styled_table(
      out, title = sprintf("Table %d. Factor Reliability (Omega)", n_table))
    .print_note(c("Total is the omega of the full scale. Omega coefficients",
                  "use the standardized solution."), width)
  }

  if (length(x$diagnostics$issues)) {
    cat("\n")
    .print_note(c("the posterior may be unreliable:",
                  paste(x$diagnostics$issues, collapse = "; ")),
                loadings$width, prefix = "Warning: ")
  }
  invisible(x)
}

# Internal functions

# This function returns the printing options of a summary, with those given
# to print() in place of those given to summary().
.loading_options <- function(x, cutoff, sort, signif_stars, name_width) {
  given <- list(cutoff = cutoff, sort = sort, signif_stars = signif_stars,
                name_width = name_width)
  out <- utils::modifyList(x$print_options, Filter(Negate(is.null), given))
  if (is.null(out$name_width)) out$name_width <- 8L
  out
}

# This function returns the variance of the standardized items that each
# factor explains, from the draws of the loadings: the posterior mean of the
# sum of their squared loadings, its share of the explained variance, and its
# share of the total variance, also accumulated. Their total equals the sum
# of the posterior means of the communalities.
.variance_table <- function(lambda, J, M) {
  ssl <- colSums(matrix(colMeans(lambda^2), J, M))
  tab <- rbind("Sum of squared loadings" = c(ssl, sum(ssl)),
               "Proportion of total" = c(ssl / sum(ssl), 1),
               "Proportion of variance" = c(ssl, sum(ssl)) / J,
               "Cumulative variance" = c(cumsum(ssl) / J, NA))
  colnames(tab) <- c(paste0("F", seq_len(M)), "Total")
  tab
}

# This function prints the tables of the loadings and the explained variance
# of a summary, numbered from `first`, with `label` after their titles. It
# returns the number of the next table and the width of the loadings.
.print_loadings <- function(x, digits, options, first, label = "") {
  df <- x$tables$table1_loadings
  df$Variable <- .short_names(df$Variable, options$name_width)
  f_cols <- colnames(df)[1L + seq_len(x$header$n_factors)]
  df <- .format_loading_table(df, x$sig_matrix, f_cols, options$cutoff,
                              options$sort, options$signif_stars, digits)$df
  cat("\n")
  width <- .print_styled_table(
    df, title = sprintf("Table %d. Standardized Factor Loadings%s", first, label))
  .print_note(c(
    sprintf("%s rotation applied; %s with %s.", x$header$rotation,
            .method_name(x$header$method), x$header$backend),
    "Diagnostics show worst-case values across factors (max Rhat, min ESS).",
    if (options$signif_stars) {
      sprintf("(*) %g%% Credible Interval excludes 0.",
              (x$probs[2L] - x$probs[1L]) * 100)
    },
    if (options$cutoff > 0) {
      sprintf("Loadings with absolute values < %.2f are hidden.", options$cutoff)
    }
  ), width)

  tab <- x$tables$table2_variance
  out <- data.frame(Statistic = rownames(tab),
                    lapply(tab, .format_numeric, digits = digits),
                    check.names = FALSE)
  cat("\n")
  table_width <- .print_styled_table(
    out, title = sprintf("Table %d. Variance Explained%s", first + 1L, label))
  .print_note(c("Posterior means. Proportion of total: share of the explained",
                "variance."), table_width)
  list(next_table = first + 2L, width = width)
}

# This function returns the fit indices of a befa_fitmeasures object as a
# matrix with their posterior mean, SD and credible interval, with Chi2_ppp
# after Chi2.
build_fit_summary_table <- function(fit_obj, probs = c(0.025, 0.975)) {
  s <- posterior::summarise_draws(
    fit_obj$fit_draws, Estimate = mean, SD = stats::sd,
    CI_Low = ~ stats::quantile(.x, probs = probs[1L], na.rm = TRUE, names = FALSE),
    CI_High = ~ stats::quantile(.x, probs = probs[2L], na.rm = TRUE, names = FALSE)
  )
  tab <- as.matrix(s[, -1L])
  rownames(tab) <- s$variable
  tab <- rbind(tab[1L, , drop = FALSE],
               Chi2_ppp = c(fit_obj$details$chi2_ppp, NA, NA, NA),
               tab[-1L, , drop = FALSE])
  colnames(tab) <- c("Estimate", "SD", "CI_Low", "CI_High")
  tab
}

# This function returns the information criteria of a befa_fitmeasures
# object with their standard error across persons and the normal interval
# that the probabilities in `probs` give.
build_information_table <- function(fit_obj, probs = c(0.025, 0.975)) {
  ic <- fit_obj$information
  tab <- cbind(ic, ic[, "Estimate"] + outer(ic[, "SE"], stats::qnorm(probs)))
  colnames(tab) <- c("Estimate", "SE", "CI_Low", "CI_High")
  tab
}

# This function prints the fit indices and the information criteria as two
# numbered tables, starting at table `first`, and returns the number of the
# next table.
.print_fit_tables <- function(fit_tab, info_tab, probs, digits, first) {
  level <- (probs[2L] - probs[1L]) * 100
  out <- as.data.frame(lapply(as.data.frame(fit_tab), .format_numeric, digits = digits))
  out <- cbind(Index = rownames(fit_tab), out)
  out[out$Index == "Chi2_ppp", c("SD", "CI_Low", "CI_High")] <- ""
  cat("\n")
  width <- .print_styled_table(
    out, title = sprintf("Table %d. Bayesian Fit Measures", first))
  .print_note(c(sprintf("Intervals are %g%% credible intervals.", level),
                "PPP: Posterior Predictive p-value (Ideal > .05)."), width)

  out <- as.data.frame(lapply(as.data.frame(info_tab), .format_numeric, digits = digits))
  out <- cbind(Criterion = rownames(info_tab), out)
  cat("\n")
  width <- .print_styled_table(
    out, title = sprintf("Table %d. Information Criteria", first + 1L))
  .print_note(c("SEs across persons (jackknife for DIC_i and p_V);",
                sprintf("%g%% normal intervals.", level)), width)
  first + 2L
}

# This function splits the omega coefficients of a befa_reliability object
# into the total score and the factors, with their posterior mean, SD and
# credible interval.
build_reliability_table <- function(rel_obj, probs = c(0.025, 0.975)) {
  s <- posterior::summarise_draws(
    rel_obj$omega_draws, Estimate = mean, SD = stats::sd,
    CI_Low = ~ stats::quantile(.x, probs = probs[1L], names = FALSE),
    CI_High = ~ stats::quantile(.x, probs = probs[2L], names = FALSE)
  )
  subscales <- as.data.frame(s[-1L, -1L])
  rownames(subscales) <- sub("^omega_", "", s$variable[-1L])
  list(total = as.data.frame(s[1L, -1L]), subscales = subscales)
}

# This function formats numbers for the tables, leaving missing values blank.
.format_numeric <- function(n, digits = 2) {
  vapply(n, function(val) {
    if (is.na(val) || identical(val, "")) return("")
    if (is.character(val)) return(val)
    sprintf(paste0("%.", digits, "f"), as.numeric(val))
  }, "")
}

# This function returns, from the draws of the loadings (with their chains),
# their posterior means (J x M), which of them have a credible interval that
# excludes zero, and their R-hat and bulk and tail ESS.
.extract_loading_summary <- function(draws_lambda, J, M, probs) {
  stacked <- .stacked_draws(draws_lambda, posterior::variables(draws_lambda))
  q_low <- apply(stacked, 2, stats::quantile, probs = probs[1L])
  q_high <- apply(stacked, 2, stats::quantile, probs = probs[2L])
  conv <- posterior::summarise_draws(draws_lambda,
                                     posterior::default_convergence_measures())
  shape <- function(x) matrix(x, nrow = J, ncol = M)
  list(lambda = shape(colMeans(stacked)),
       sig_matrix = shape((q_low > 0) | (q_high < 0)),
       rhat = shape(conv$rhat), ess_bulk = shape(conv$ess_bulk),
       ess_tail = shape(conv$ess_tail))
}

# This function formats the table of loadings: sorted by primary factor if
# asked, with significance stars and hidden small loadings.
.format_loading_table <- function(df, sig_mat, f_cols, cutoff, sort,
                                  signif_stars, digits) {
  if (sort) {
    loadings <- abs(as.matrix(df[, f_cols, drop = FALSE]))
    ord <- order(apply(loadings, 1, which.max), -apply(loadings, 1, max))
    df <- df[ord, ]
    sig_mat <- sig_mat[ord, , drop = FALSE]
  }
  for (col in setdiff(names(df), "Variable")) {
    vals <- df[[col]]
    if (col %in% f_cols) {
      formatted <- .format_numeric(vals, digits)
      if (signif_stars) {
        star <- sig_mat[, match(col, f_cols)]
        formatted <- paste0(formatted, ifelse(star, "*", " "))
      }
      if (cutoff > 0) formatted[abs(as.numeric(vals)) < cutoff] <- ""
      df[[col]] <- formatted
    } else if (col %in% c("EssBulk", "EssTail")) {
      df[[col]] <- sprintf("%.0f", vals)
    } else {
      df[[col]] <- .format_numeric(vals, digits)
    }
  }
  list(df = df, sig_mat = sig_mat)
}

# This function prints a table between thick lines, with the first column on
# the left and the others on the right, and returns its width.
.print_styled_table <- function(df, title = NULL) {
  cn <- colnames(df)
  widths <- vapply(seq_along(df), function(i) {
    max(nchar(cn[i]), max(nchar(df[[i]])))
  }, 1L) + 2L
  total_width <- sum(widths)
  thick_line <- paste0(rep("\u2017", total_width), collapse = "")
  thin_line <- paste0(rep("\u2014", total_width), collapse = "")
  cell <- function(text, j) {
    format(text, width = widths[j], justify = if (j == 1L) "left" else "right")
  }

  if (!is.null(title)) cat(title, "\n")
  cat(thick_line, "\n")
  for (j in seq_along(df)) cat(cell(cn[j], j))
  cat("\n", thin_line, "\n", sep = "")
  for (i in seq_len(nrow(df))) {
    for (j in seq_along(df)) cat(cell(df[i, j], j))
    cat("\n")
  }
  cat(thick_line, "\n")
  invisible(total_width)
}

# This function prints the note under a table, wrapped to its width.
.print_note <- function(notes, width, prefix = "Note: ") {
  text <- paste0(prefix, paste(notes, collapse = " "))
  cat(paste(strwrap(text, width = width), collapse = "\n"), "\n")
}

# This function shortens the names longer than `width` characters by
# cutting their middle, which keeps their start and their end, where item
# numbers usually are. Names that would become equal end with their item
# number instead.
.short_names <- function(x, width = 8L) {
  long <- nchar(x) > width
  if (!any(long)) return(x)
  head <- ceiling((width - 1L) / 2)
  tail <- width - 1L - head
  out <- x
  out[long] <- paste0(substr(x[long], 1L, head), "\u2026",
                      substring(x[long], nchar(x[long]) - tail + 1L))
  clash <- long & (duplicated(out) | duplicated(out, fromLast = TRUE))
  index <- as.character(which(clash))
  out[clash] <- paste0(substr(x[clash], 1L, width - 1L - nchar(index)), "\u2026",
                       index)
  out
}
