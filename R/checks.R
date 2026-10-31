# checks.R
# Checks of the arguments that the internal functions receive
# Author: Ricardo Rey-Sáez
# Last modified: 31-10-2026

# This function checks that a value is one whole number of at least `lower`
# (a positive one by default).
.check_count <- function(x, arg, lower = 1) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < lower ||
      x != round(x)) {
    .befa_stop("`", arg, "` must be ",
               if (lower == 1) "a positive" else "a non-negative",
               " whole number.")
  }
  invisible(TRUE)
}

# This function checks that a value is one number strictly between 0 and 1.
.check_probability <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= 0 || x >= 1) {
    .befa_stop("`", arg, "` must be a number between 0 and 1.")
  }
  invisible(TRUE)
}

# This function checks that a value is one of a few strings and returns it.
.check_choice <- function(x, choices, arg) {
  if (!is.character(x) || length(x) != 1L || !x %in% choices) {
    .befa_stop("`", arg, "` must be one of ",
               paste0("\"", choices, "\"", collapse = ", "), ".")
  }
  x
}

# This function checks the arguments of befa() that the later steps do not
# check themselves: the logical flags, one source of data (raw data or
# summary statistics), FIML only with raw data, and numeric data.
.check_befa_args <- function(data, n_factors, ordered, sample_nobs, sample_cov,
                             sample_cor, missing, factor_scores,
                             compute_fit_indices, compute_reliability,
                             verbose) {
  flags <- list(ordered = ordered, factor_scores = factor_scores,
                compute_fit_indices = compute_fit_indices,
                compute_reliability = compute_reliability, verbose = verbose)
  for (arg in names(flags)) .check_flag(flags[[arg]], arg)
  if (ordered) {
    .befa_stop("Ordinal factor analysis (`ordered = TRUE`) is not implemented ",
               "yet; use `ordered = FALSE`.")
  }
  if (is.null(n_factors)) .befa_stop("`n_factors` is required.")
  .check_count(n_factors, "n_factors")

  has_data <- !is.null(data)
  has_summary <- !is.null(sample_nobs) || !is.null(sample_cov) ||
    !is.null(sample_cor)
  if (has_data && has_summary) {
    .befa_stop("Provide either `data` or summary statistics, not both.")
  }
  if (!has_data && !has_summary) {
    .befa_stop("Provide `data`, or summary statistics (`sample_nobs` with ",
               "`sample_cov` or `sample_cor`).")
  }
  if (missing == "FIML" && !has_data) {
    .befa_stop("FIML needs the raw `data`, not summary statistics.")
  }
  if (has_data) {
    if (!is.data.frame(data) && !is.matrix(data)) {
      .befa_stop("`data` must be a data frame or a matrix.")
    }
    if (!all(vapply(as.data.frame(data), is.numeric, TRUE))) {
      .befa_stop("`data` must contain numeric columns only.")
    }
  } else {
    if (!is.null(sample_cov) && !is.null(sample_cor)) {
      .befa_stop("Provide either `sample_cov` or `sample_cor`, not both.")
    }
    for (arg in c("sample_cov", "sample_cor")) {
      S <- list(sample_cov = sample_cov, sample_cor = sample_cor)[[arg]]
      if (is.null(S)) next
      .check_matrix(S, arg)
      if (nrow(S) != ncol(S) || !isSymmetric(unname(S))) {
        .befa_stop("`", arg, "` must be a symmetric square matrix.")
      }
    }
  }
  invisible(TRUE)
}

# This function checks that a value is TRUE or FALSE.
.check_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    .befa_stop("`", arg, "` must be TRUE or FALSE.")
  }
  invisible(TRUE)
}

# This function checks that a value is a numeric matrix.
.check_matrix <- function(x, arg) {
  if (!is.matrix(x) || !is.numeric(x)) {
    .befa_stop("`", arg, "` must be a numeric matrix.")
  }
  invisible(TRUE)
}
