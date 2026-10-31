# checks.R
# Checks of the arguments that the internal functions receive
# Author: Ricardo Rey-Sáez
# Last modified: 31-10-2026

# This function checks that a value is one positive whole number.
.check_count <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 1 ||
      x != round(x)) {
    .befa_stop("`", arg, "` must be a positive whole number.")
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
