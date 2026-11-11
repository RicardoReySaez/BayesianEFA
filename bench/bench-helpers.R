# bench-helpers.R
# Shared helpers of the benchmarks (not part of the package build)
# Author: Ricardo Rey-Sáez
# Last modified: 31-10-2026

# The benchmarks load the package from source to reach its internal
# functions, and are run from the package root, e.g.
#   Rscript bench/bench-engines.R

suppressMessages(pkgload::load_all(".", quiet = TRUE))

# This function simulates raw data from an M-factor model with simple
# structure: item j loads 0.65 on factor ((j - 1) %% M) + 1 and 0.1 on the
# others, the uniquenesses complete unit variances, and a share `missing` of
# the values is removed completely at random (no row is left empty).
.bench_data <- function(J, M, N, missing = 0, seed = 1) {
  set.seed(seed)
  L <- matrix(0.1, J, M)
  L[cbind(seq_len(J), (seq_len(J) - 1L) %% M + 1L)] <- 0.65
  L <- L / pmax(1, sqrt(rowSums(L^2)) / 0.9)
  Sigma <- tcrossprod(L) + diag(1 - rowSums(L^2))
  Y <- mvnfast::rmvn(N, mu = rep(0, J), sigma = Sigma)
  if (missing > 0) {
    Y[sample(N * J, round(missing * N * J))] <- NA
    empty <- rowSums(!is.na(Y)) == 0L
    Y[empty, 1L] <- 0
  }
  Y
}

# This function returns what one benchmark condition needs to run Stan: the
# structure, its code, the data passed to Stan and one list of starting
# values, for the correlation model with the unit-vector prior.
.bench_model <- function(Y, M, engine, model = "cor") {
  missing <- if (anyNA(Y)) "FIML" else "listwise"
  prepared <- .prepare_befa_data(Y, model, missing)
  spec <- .befa_spec(model, "unit_vector", M, .data_patterns(prepared),
                     engine = engine)
  prior <- prepare_befa_priors(list(), model, "unit_vector", M, prepared$J)
  list(spec = spec, code = .stan_code(spec),
       data = .stan_data(spec, prepared, M, prior),
       init = .start_values(spec, prepared, M, prior, 1L, 1L)[[1L]],
       patterns = length(prepared$groups))
}

# This function appends rows to a CSV file of results, writing its header
# the first time.
.bench_write <- function(row, file) {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  utils::write.table(row, file, sep = ",", row.names = FALSE,
                     col.names = !file.exists(file), append = file.exists(file))
}
