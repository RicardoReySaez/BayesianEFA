# bench-backends.R
# The same models with each backend and method, to set backend = "auto"
# Author: Ricardo Rey-Sáez
# Last modified: 31-10-2026

# Each backend starts from an empty cache, so the first fit includes the
# compilation (none for stanli) and the second one reuses the compiled
# program. The fits use the defaults of befa(). ADVI and Pathfinder are then
# compared with MCMC on the communalities.

source("bench/bench-helpers.R")
out <- "bench/results/backends.csv"
out_methods <- "bench/results/methods.csv"
unlink(c(out, out_methods))

# This function fits a model and returns it with its wall time in seconds.
.timed_befa <- function(...) {
  start <- Sys.time()
  fit <- suppressWarnings(befa(..., verbose = FALSE, compute_fit_indices = FALSE,
                               compute_reliability = FALSE))
  list(fit = fit, wall = as.numeric(difftime(Sys.time(), start, units = "secs")))
}

datasets <- list(
  HS_M3 = list(Y = HS_data, M = 3L),
  sim_J30_M3_N500 = list(Y = .bench_data(30, 3, 500, seed = 11), M = 3L)
)
mcmc <- list()
for (backend in c("cmdstanr", "stanli", "stanr", "rstan")) {
  options(BayesianEFA.cache_dir = file.path(tempdir(), paste0("cache-", backend)))
  for (name in names(datasets)) {
    d <- datasets[[name]]
    first <- .timed_befa(d$Y, d$M, backend = backend, seed = 1)
    second <- .timed_befa(d$Y, d$M, backend = backend, seed = 1)
    fit <- second$fit
    ess <- posterior::summarise_draws(
      posterior::subset_draws(fit$draws, variable = c("h2", "Psi_std")),
      "ess_bulk")$ess_bulk
    row <- data.frame(
      backend = backend, data = name, engine = fit$spec$engine,
      wall_first = first$wall, wall_cached = second$wall,
      chain_seconds = sum(fit$fit$seconds), min_ess_bulk = min(ess),
      ess_per_chain_second = min(ess) / sum(fit$fit$seconds),
      divergent = fit$diagnostics$divergent,
      max_rhat = fit$diagnostics$max_rhat
    )
    .bench_write(row, out)
    cat(sprintf("%-8s %-16s first %6.1f s  cached %5.1f s  ESS/s %7.1f\n",
                backend, name, first$wall, second$wall,
                row$ess_per_chain_second))
    mcmc[[paste(backend, name)]] <- fit
  }
}

# Approximations against the cmdstanr MCMC fit of each data set
for (backend in c("cmdstanr", "stanli")) {
  options(BayesianEFA.cache_dir = file.path(tempdir(), paste0("cache-", backend)))
  for (name in names(datasets)) {
    d <- datasets[[name]]
    ref <- colMeans(.stacked_draws(mcmc[[paste("cmdstanr", name)]]$draws, "h2"))
    for (method in c("advi", "pathfinder")) {
      res <- .timed_befa(d$Y, d$M, backend = backend, method = method, seed = 1)
      h2 <- colMeans(.stacked_draws(res$fit$draws, "h2"))
      row <- data.frame(backend = backend, data = name, method = method,
                        wall = res$wall, max_abs_error_h2 = max(abs(h2 - ref)),
                        pareto_k = res$fit$diagnostics$pareto_k)
      .bench_write(row, out_methods)
      cat(sprintf("%-8s %-16s %-10s %5.1f s  max|dh2| %.3f  k %.2f\n", backend,
                  name, method, res$wall, row$max_abs_error_h2, row$pareto_k))
    }
  }
}
