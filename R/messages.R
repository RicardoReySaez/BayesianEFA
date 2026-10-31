# messages.R
# Errors, warnings and messages shown by the package
# Author: Ricardo Rey-Sáez
# Last modified: 31-10-2026

# This function stops the analysis with an error.
.befa_stop <- function(...) {
  stop(paste0(...), call. = FALSE)
}

# This function shows a warning.
.befa_warn <- function(...) {
  warning(paste0(...), call. = FALSE)
}

# This function explains a choice made by the package.
.befa_msg <- function(...) {
  message(paste0("BayesianEFA: ", ...))
}
