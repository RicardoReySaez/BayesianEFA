# bench-backends-large.R
# cmdstanr and stanli on a larger model
# Author: Ricardo Rey-Sáez
# Last modified: 01-11-2026

# stanli interprets the program as a graph of Stan Math operations, so its
# cost per operation is higher than that of compiled code. This fit checks
# whether that matters on a larger model. The compiled program comes from the
# cache of the other benchmarks, so the times exclude compilation. (FIML with
# many patterns is compared per gradient in bench-engines.R; full fits of it
# take too long to be worth running.)

source("bench/bench-helpers.R")
out <- "bench/results/backends-large.csv"
unlink(out)
options(BayesianEFA.cache_dir = file.path("bench", ".exe", "cache"))

datasets <- list(
  sim_J60_M4_N1000 = list(Y = .bench_data(60, 4, 1000, seed = 21), M = 4L)
)
for (name in names(datasets)) {
  d <- datasets[[name]]
  missing <- if (is.null(d$missing)) "listwise" else d$missing
  for (backend in c("cmdstanr", "stanli")) {
    # A first fit compiles the program if the cache does not have it yet.
    if (backend == "cmdstanr") {
      suppressWarnings(befa(d$Y, d$M, missing = missing, backend = backend,
                            chains = 1, iter_warmup = 10, iter_sampling = 10,
                            verbose = FALSE, compute_fit_indices = FALSE,
                            compute_reliability = FALSE))
    }
    start <- Sys.time()
    fit <- suppressWarnings(befa(d$Y, d$M, missing = missing, backend = backend,
                                 seed = 1, verbose = FALSE,
                                 compute_fit_indices = FALSE,
                                 compute_reliability = FALSE))
    wall <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    ess <- posterior::summarise_draws(
      posterior::subset_draws(fit$draws, variable = c("h2", "Psi_std")),
      "ess_bulk")$ess_bulk
    leapfrog <- sum(.stacked_draws(fit$fit$sampler, "n_leapfrog__"))
    row <- data.frame(
      backend = backend, data = name, engine = fit$spec$engine, wall = wall,
      chain_seconds = sum(fit$fit$seconds), min_ess_bulk = min(ess),
      ess_per_chain_second = min(ess) / sum(fit$fit$seconds),
      us_per_leapfrog = sum(fit$fit$seconds) / leapfrog * 1e6
    )
    .bench_write(row, out)
    cat(sprintf("%-8s %-24s wall %6.1f s  ESS/s %6.1f  %7.1f us/leapfrog\n",
                backend, name, wall, row$ess_per_chain_second,
                row$us_per_leapfrog))
  }
}
