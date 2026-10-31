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

# This function checks that a value is a numeric matrix.
.check_matrix <- function(x, arg) {
  if (!is.matrix(x) || !is.numeric(x)) {
    .befa_stop("`", arg, "` must be a numeric matrix.")
  }
  invisible(TRUE)
}
