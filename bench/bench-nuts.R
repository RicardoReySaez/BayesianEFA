# bench-nuts.R
# NUTS with each target acceptance rate, to review the defaults of befa()
# Author: Ricardo Rey-Sáez
# Last modified: 31-10-2026

# Fits three data sets with cmdstanr and adapt_delta in {0.8, 0.9, 0.95},
# twice each, and records the divergent transitions, the tree depths, the
# leapfrog steps per iteration, the smallest bulk ESS of h2 and Psi_std (which
# do not depend on the rotation) and that ESS per second of sampling, summed
# over chains.

source("bench/bench-helpers.R")
out <- "bench/results/nuts.csv"
unlink(out)
options(BayesianEFA.cache_dir = file.path("bench", ".exe", "cache"))

datasets <- list(
  HS_M3 = list(Y = HS_data, M = 3L),
  sim_J30_M3_N500 = list(Y = .bench_data(30, 3, 500, seed = 11), M = 3L),
  sim_J12_M4_N150 = list(Y = .bench_data(12, 4, 150, seed = 12), M = 4L)
)
for (name in names(datasets)) {
  d <- datasets[[name]]
  for (delta in c(0.8, 0.9, 0.95)) {
    for (seed in 1:2) {
      fit <- suppressWarnings(befa(
        d$Y, d$M, backend = "cmdstanr", adapt_delta = delta, seed = seed,
        verbose = FALSE, compute_fit_indices = FALSE,
        compute_reliability = FALSE))
      s <- fit$fit$sampler
      depth <- .stacked_draws(s, "treedepth__")
      leapfrog <- .stacked_draws(s, "n_leapfrog__")
      ess <- posterior::summarise_draws(
        posterior::subset_draws(fit$draws, variable = c("h2", "Psi_std")),
        "ess_bulk")$ess_bulk
      seconds <- sum(fit$fit$seconds)
      row <- data.frame(
        data = name, adapt_delta = delta, seed = seed,
        engine = fit$spec$engine, divergent = fit$diagnostics$divergent,
        max_treedepth = max(depth), share_depth_10 = mean(depth >= 10),
        mean_leapfrog = mean(leapfrog),
        stepsize = mean(.stacked_draws(s, "stepsize__")),
        chain_seconds = seconds, min_ess_bulk = min(ess),
        ess_per_second = min(ess) / seconds,
        max_rhat = fit$diagnostics$max_rhat
      )
      .bench_write(row, out)
      cat(sprintf("%-16s delta=%.2f seed=%d div=%3d depth=%2d leap=%5.1f ESS/s=%7.1f\n",
                  name, delta, seed, row$divergent, row$max_treedepth,
                  row$mean_leapfrog, row$ess_per_second))
    }
  }
}
