# bench-engines.R
# Cost of one gradient of the log density with each likelihood engine
# Author: Ricardo Rey-Sáez
# Last modified: 31-10-2026

# The programs are compiled once with cmdstanr. Each condition runs NUTS
# without adaptation, with a tiny step size and a maximum tree depth of 5, so
# every iteration takes the same 31 leapfrog steps; the sampling time that
# CmdStan reports, divided by the leapfrog steps, is the cost of one gradient
# (process start and data reading excluded). The number of iterations is set
# from a short probe so that each condition samples for about 1.5 seconds.

source("bench/bench-helpers.R")
out <- "bench/results/engines.csv"
unlink(out)
exe_dir <- file.path("bench", ".exe")
dir.create(exe_dir, showWarnings = FALSE)

# This function compiles (once) the program of a structure with cmdstanr.
.bench_compile <- local({
  cache <- list()
  function(spec, code) {
    label <- .spec_label(spec)
    if (is.null(cache[[label]])) {
      file <- cmdstanr::write_stan_file(code, dir = exe_dir,
                                        basename = gsub("-", "_", label))
      cache[[label]] <<- cmdstanr::cmdstan_model(file, quiet = TRUE)
    }
    cache[[label]]
  }
})

# This function returns the seconds per gradient of `iter` iterations.
.bench_run <- function(model, bm, iter) {
  fit <- suppressMessages(suppressWarnings(model$sample(
    data = bm$data, chains = 1, iter_warmup = 0, iter_sampling = iter,
    adapt_engaged = FALSE, step_size = 1e-4, max_treedepth = 5,
    init = list(bm$init), seed = 1, refresh = 0, show_messages = FALSE,
    show_exceptions = FALSE)))
  leapfrog <- sum(fit$sampler_diagnostics(format = "draws_matrix")[, "n_leapfrog__"])
  fit$time()$chains$sampling / leapfrog
}

grid <- expand.grid(J = c(10, 15, 20, 30, 60, 100), M = c(2, 3, 5, 8),
                    missing = c(0, 0.1, 0.3), KEEP.OUT.ATTRS = FALSE)
grid <- grid[(grid$J - grid$M)^2 >= grid$J + grid$M, ]
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  Y <- .bench_data(g$J, g$M, N = 1000, missing = g$missing, seed = i)
  for (engine in c("cholesky", "woodbury")) {
    bm <- .bench_model(Y, g$M, engine)
    model <- .bench_compile(bm$spec, bm$code)
    probe <- .bench_run(model, bm, 2)
    iter <- max(2, min(300, round(1.5 / (probe * 31))))
    sec <- .bench_run(model, bm, iter)
    row <- data.frame(J = g$J, M = g$M, N = 1000, missing = g$missing,
                      patterns = bm$patterns, engine = engine, iter = iter,
                      us_per_gradient = sec * 1e6)
    .bench_write(row, out)
    cat(sprintf("J=%3d M=%d missing=%.1f P=%4d %-8s %10.1f us\n", g$J, g$M,
                g$missing, bm$patterns, engine, sec * 1e6))
  }
}
