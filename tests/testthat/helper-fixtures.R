# helper-fixtures.R
# Shared helpers for the tests that compare with BayesianEFA v0.1.0
# Author: Ricardo Rey-Sáez
# Last modified: 31-10-2026

# This function reads a golden fixture that notes/fixtures_v010.R made from
# the v0.1.0 sources: a list with the inputs, the output and its metadata.
.read_fixture <- function(name) {
  readRDS(testthat::test_path("fixtures", paste0(name, "_v010.rds")))
}

# This function rebuilds the full data matrix from the complete and incomplete
# rows that the FIML helpers of v0.1.0 take.
.fiml_data_matrix <- function(data_fiml) {
  J <- ncol(data_fiml$Y_comp)
  Y <- matrix(NA_real_, length(data_fiml$comp_idx) + length(data_fiml$miss_idx), J)
  Y[data_fiml$comp_idx, ] <- data_fiml$Y_comp
  Y[data_fiml$miss_idx, ] <- data_fiml$Y_miss
  Y
}
