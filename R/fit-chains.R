# fit-chains.R
# Run Stan in background processes behind one progress display
# Author: Ricardo Rey-Sáez
# Last modified: 31-10-2026

# Internal functions

# This function fits a program with one backend and method, and returns its
# draws in a common format (class `befa_draws`). The program is compiled once
# and kept in the cache. Each MCMC chain runs in its own background R process,
# at most `parallel_chains` at a time, with its own starting values and seed;
# ADVI and Pathfinder run in one process. `inits` holds one list of starting
# values per chain (or per Pathfinder path).
.run_stan <- function(code, data, inits, backend, settings, quiet = FALSE) {
  files <- .compiled_model(code, backend, quiet = quiet)
  show <- !quiet && settings$refresh > 0
  if (show) .stan_header(backend, settings)

  job <- function(init, chain) {
    list(func = .stan_job,
         args = list(backend = backend, method = settings$method,
                     stan_file = files$stan_file, data = data,
                     settings = settings, init = init, chain = chain))
  }
  if (settings$method == "mcmc") {
    chains <- settings$chains
    jobs <- lapply(seq_len(chains), function(ch) job(inits[[ch]], ch))
    label <- if (chains == 1L) "1 chain" else paste(chains, "chains")
    results <- .run_jobs(jobs, cores = settings$parallel_chains,
                         per_job = settings$iter_warmup + settings$iter_sampling,
                         warmup = settings$iter_warmup, label = "Chain",
                         done = paste("Sampled", label), quiet = !show)
  } else {
    results <- .run_jobs(list(job(inits, 1L)), label = .method_name(settings$method),
                         done = paste("Ran", .method_name(settings$method)),
                         quiet = !show)
  }
  fit <- .combine_jobs(results, backend, settings)
  if (settings$method == "pathfinder" && backend %in% c("cmdstanr", "stanr")) {
    fit <- .pathfinder_resample(fit, settings)
  }
  fit
}

# This function finishes multi-path Pathfinder: it records the Pareto k of
# the importance ratios of all the candidate draws and, with
# `psis_resample = TRUE`, resamples `draws` of them with their PSIS weights.
# CmdStan could resample itself, but its resampled draws repeat, and the
# ties in the tail leave the Pareto k undefined.
.pathfinder_resample <- function(fit, settings) {
  if (!all(c("lp__", "lp_approx__") %in% posterior::variables(fit$draws))) {
    return(fit)
  }
  ratios <- as.vector(.stacked_draws(fit$draws, "lp__") -
                        .stacked_draws(fit$draws, "lp_approx__"))
  if (!all(is.finite(ratios))) return(fit)
  psis <- suppressWarnings(loo::psis(ratios, r_eff = 1))
  fit$pareto_k <- as.numeric(psis$diagnostics$pareto_k)
  if (isTRUE(settings$psis_resample)) {
    weights <- exp(as.vector(stats::weights(psis, log = TRUE)))
    fit$draws <- withr::with_seed(settings$seed, posterior::resample_draws(
      posterior::merge_chains(fit$draws), weights = weights,
      method = "stratified", ndraws = settings$draws))
  }
  fit
}

# This function joins the results of the jobs into a `befa_draws` object: the
# draws of the parameters with lp__ (and lp_approx__ for the approximations),
# the sampler diagnostics of MCMC, and the time each job took.
.combine_jobs <- function(results, backend, settings) {
  draws <- lapply(results, function(r) .standardize_draws(r$draws))
  sampler <- lapply(results, `[[`, "sampler")
  structure(list(
    draws = if (length(draws) == 1L) draws[[1L]] else posterior::bind_draws(draws, along = "chain"),
    sampler = if (all(lengths(sampler))) posterior::bind_draws(sampler, along = "chain"),
    method = settings$method, backend = backend, settings = settings,
    seconds = vapply(results, `[[`, 0, "seconds")
  ), class = "befa_draws")
}

# This function gives the extra variables of every backend the same names:
# lp__ for the log density of the target and lp_approx__ for that of an
# approximation (rstan writes them as log_p__ and log_g__). Other bookkeeping
# variables, such as Pathfinder's path__, are dropped.
.standardize_draws <- function(draws) {
  draws <- posterior::as_draws_array(draws)
  vars <- posterior::variables(draws)
  rename <- c(log_p__ = "lp__", log_g__ = "lp_approx__")
  if ("log_p__" %in% vars) {
    draws <- posterior::subset_draws(draws, variable = setdiff(vars, "lp__"))
    vars <- posterior::variables(draws)
  }
  hit <- vars %in% names(rename)
  vars[hit] <- rename[vars[hit]]
  posterior::variables(draws) <- vars
  keep <- !grepl("__$", vars) | vars %in% c("lp__", "lp_approx__")
  posterior::subset_draws(draws, variable = vars[keep])
}

# This function announces the fit: the backend and how the chains run.
.stan_header <- function(backend, settings) {
  cli::cli_rule(left = "Bayesian exploratory factor analysis", right = backend)
  s <- settings
  text <- switch(s$method,
    mcmc = paste0(s$chains, " chain", if (s$chains > 1L) "s", " of ",
                  s$iter_sampling, " draws after ", s$iter_warmup,
                  " warmup iterations, ", min(s$parallel_chains, s$chains),
                  " at a time."),
    advi = paste0("ADVI (", s$algorithm, ") with up to ", s$iter,
                  " iterations and ", s$draws, " draws."),
    pathfinder = {
      # stanli runs a single path
      paths <- if (backend == "stanli") 1L else s$num_paths
      paste0("Pathfinder with ", paths, " path", if (paths > 1L) "s",
             " and ", s$draws, " draws.")
    }
  )
  cli::cli_alert_info(text)
}

# This function fits one chain (or one ADVI or Pathfinder run) in the
# background process that runs it, and returns its draws as posterior
# objects. It must be self-contained: the process does not load BayesianEFA,
# only the backend and posterior. Each backend reuses the model compiled in
# the folder of `stan_file`.
.stan_job <- function(backend, method, stan_file, data, settings, init, chain) {
  s <- settings
  start <- Sys.time()
  as_array <- function(x) posterior::as_draws_array(x)
  sampler <- NULL

  if (backend == "rstan") {
    model <- rstan::stan_model(file = stan_file, auto_write = TRUE)
    if (method == "mcmc") {
      fit <- rstan::sampling(model, data = data, chains = 1L, chain_id = chain,
                             iter = s$iter_warmup + s$iter_sampling,
                             warmup = s$iter_warmup, seed = s$seed,
                             init = list(init), refresh = s$refresh,
                             control = list(adapt_delta = s$adapt_delta,
                                            max_treedepth = s$max_treedepth))
      sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)[[1L]]
      sampler <- as_array(array(sp, dim = c(nrow(sp), 1L, ncol(sp)),
                                dimnames = list(NULL, NULL, colnames(sp))))
    } else {
      fit <- rstan::vb(model, data = data, algorithm = s$algorithm,
                       iter = s$iter, output_samples = s$draws,
                       tol_rel_obj = s$tol_rel_obj, seed = s$seed, init = init,
                       refresh = s$refresh)
    }
    draws <- as_array(fit)
  } else if (backend == "cmdstanr") {
    model <- cmdstanr::cmdstan_model(stan_file, quiet = TRUE)
    fit <- switch(method,
      mcmc = model$sample(data = data, chains = 1L, chain_ids = chain,
                          iter_warmup = s$iter_warmup,
                          iter_sampling = s$iter_sampling, seed = s$seed,
                          init = list(init), refresh = s$refresh,
                          adapt_delta = s$adapt_delta,
                          max_treedepth = s$max_treedepth),
      advi = model$variational(data = data, algorithm = s$algorithm,
                               iter = s$iter, draws = s$draws,
                               tol_rel_obj = s$tol_rel_obj, seed = s$seed,
                               init = list(init), refresh = s$refresh),
      # Every candidate draw comes back; .pathfinder_resample() resamples them.
      pathfinder = model$pathfinder(data = data, draws = s$draws,
                                    num_paths = s$num_paths,
                                    max_lbfgs_iters = s$iter,
                                    psis_resample = FALSE,
                                    seed = s$seed, init = init,
                                    refresh = s$refresh)
    )
    draws <- as_array(fit$draws())
    if (method == "mcmc") sampler <- as_array(fit$sampler_diagnostics())
  } else if (backend == "stanr") {
    model <- stanr::stan_model(stan_file = stan_file, quiet = TRUE)
    if (method == "pathfinder") init <- init[[1L]]
    fit <- switch(method,
      mcmc = model$sample(data = data, chains = 1L, chain_ids = chain,
                          iter_warmup = s$iter_warmup,
                          iter_sampling = s$iter_sampling, seed = s$seed,
                          init = init, refresh = s$refresh,
                          adapt_delta = s$adapt_delta,
                          max_treedepth = s$max_treedepth),
      advi = model$variational(data = data, algorithm = s$algorithm,
                               iter = s$iter, draws = s$draws,
                               tol_rel_obj = s$tol_rel_obj, seed = s$seed,
                               init = init, refresh = s$refresh),
      # Every candidate draw comes back; .pathfinder_resample() resamples them.
      pathfinder = model$pathfinder(data = data, draws = s$draws,
                                    num_paths = s$num_paths,
                                    max_lbfgs_iters = s$iter,
                                    psis_resample = FALSE,
                                    seed = s$seed, init = init,
                                    refresh = s$refresh)
    )
    draws <- as_array(fit$draws())
    if (method == "mcmc") sampler <- as_array(fit$sampler_diagnostics())
  } else {
    model <- stanli::stanli_model(file = stan_file, data = data, seed = s$seed)
    if (method == "pathfinder") init <- init[[1L]]
    q <- stanli::unconstrain(model, init)
    fit <- switch(method,
      mcmc = stanli::sample_model(model, chains = 1L, chain_id = chain,
                                  warmup = s$iter_warmup,
                                  samples = s$iter_sampling, seed = s$seed,
                                  init = q, refresh = s$refresh,
                                  delta = s$adapt_delta,
                                  max_depth = s$max_treedepth),
      advi = stanli::variational_model(model, algorithm = s$algorithm,
                                       iter = s$iter, draws = s$draws,
                                       tol_rel_obj = s$tol_rel_obj,
                                       seed = s$seed, init = q,
                                       refresh = s$refresh),
      # stanli runs a single Pathfinder path, without PSIS resampling
      pathfinder = stanli::pathfinder_model(model, draws = s$draws,
                                            max_lbfgs_iters = s$iter,
                                            seed = s$seed, init = q,
                                            refresh = s$refresh)
    )
    all <- stanli::as_draws_array(fit, include_sampler = TRUE)
    vars <- posterior::variables(all)
    diag <- c("accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__",
              "divergent__", "energy__")
    draws <- posterior::subset_draws(all, variable = setdiff(vars, diag))
    if (method == "mcmc") {
      sampler <- posterior::subset_draws(all, variable = intersect(vars, diag))
    }
  }

  list(draws = draws, sampler = sampler,
       seconds = as.numeric(difftime(Sys.time(), start, units = "secs")))
}

# This function runs jobs in background R processes, at most `cores` at a
# time, and returns what each returned. Stan's messages stay in the files of
# those processes; with `per_job` (the iterations of each MCMC chain) the
# iterations they report move a single progress bar, and otherwise a spinner
# shows the time elapsed. The bar redraws itself in place only in a console
# that can do that; in knitr documents and logs a single line reports the time
# taken. It starts with the first progress report, so the seconds the
# processes take to start do not distort its estimate of the time left.
# Interrupting the fit stops every process.
.run_jobs <- function(jobs, cores = length(jobs), per_job = NULL, warmup = 0L,
                      label = "Job", done = "Finished", quiet = FALSE) {
  n <- length(jobs)
  results <- vector("list", n)
  progress <- integer(n)
  queue <- seq_len(n)
  running <- list()
  on.exit(for (job in running) .stop_job(job), add = TRUE)

  dynamic <- !quiet && cli::is_dynamic_tty()
  started <- Sys.time()
  bar <- NULL
  if (dynamic) {
    old <- options(cli.progress_show_after = 0,
                   cli.progress_bar_style = .bar_style())
    on.exit(options(old), add = TRUE)
    if (is.null(per_job)) {
      bar <- cli::cli_progress_bar(
        total = NA, clear = FALSE, auto_terminate = FALSE,
        format = paste("{cli::pb_spin}", label, "| {cli::pb_elapsed}"),
        format_done = paste("{cli::col_green(cli::symbol$tick)}", done,
                            "in {cli::pb_elapsed}."))
    }
  }

  while (length(queue) || length(running)) {
    while (length(queue) && length(running) < cores) {
      i <- queue[1L]
      running[[as.character(i)]] <- .launch_job(jobs[[i]]$func, jobs[[i]]$args)
      queue <- queue[-1L]
    }
    for (key in names(running)) {
      job <- running[[key]]
      i <- as.integer(key)
      alive <- job$process$is_alive()
      if (!is.null(per_job)) progress[i] <- max(progress[i], .reported_iteration(job))
      if (!alive) {
        results[i] <- list(.job_value(job, if (n > 1L) paste(label, i) else label))
        if (!is.null(per_job)) progress[i] <- per_job
        .stop_job(job)
        running[[key]] <- NULL
      }
    }
    if (dynamic && !is.null(per_job)) {
      phase <- if (all(progress > warmup)) "Sampling" else "Warmup  "
      if (is.null(bar) && max(progress) > 1L) {
        bar <- cli::cli_progress_bar(
          total = n * per_job, clear = FALSE, auto_terminate = FALSE,
          status = phase,
          format = paste("{cli::pb_spin} {cli::pb_status} {cli::pb_bar}",
                         "{cli::pb_percent} | ETA: {cli::pb_eta}"),
          format_done = paste("{cli::col_green(cli::symbol$tick)}", done,
                              "in {cli::pb_elapsed}."))
      }
      if (!is.null(bar)) {
        cli::cli_progress_update(id = bar, set = sum(progress), status = phase)
      }
    } else if (!is.null(bar)) {
      cli::cli_progress_update(id = bar)
    }
    if (length(running)) Sys.sleep(0.1)
  }

  if (!is.null(bar)) {
    cli::cli_progress_done(id = bar)
  } else if (!quiet) {
    elapsed <- .format_elapsed(as.numeric(difftime(Sys.time(), started,
                                                   units = "secs")))
    cli::cli_alert_success("{done} in {elapsed}.")
  }
  results
}

# This function starts `func(args)` in a background Rscript process, which
# saves its value (or its error) to a file. The process writes Stan's output
# to files that the user's session reads for the progress bar, and loads only
# what the job needs, so no backend is ever loaded in the user's session.
# `func` must be self-contained, as the process does not load BayesianEFA.
.launch_job <- function(func, args) {
  environment(func) <- globalenv()
  files <- list(job = tempfile("befa-job-", fileext = ".rds"),
                result = tempfile("befa-result-", fileext = ".rds"),
                out = tempfile("befa-out-", fileext = ".txt"),
                err = tempfile("befa-err-", fileext = ".txt"),
                script = tempfile("befa-run-", fileext = ".R"))
  saveRDS(list(func = func, args = args), files$job)
  writeLines(c(
    "a <- commandArgs(TRUE)",
    "job <- readRDS(a[1L])",
    "res <- tryCatch(list(ok = TRUE, value = do.call(job$func, job$args)),",
    "                error = function(e) list(ok = FALSE, message = conditionMessage(e)))",
    "saveRDS(res, a[2L])"
  ), files$script)
  process <- callr::rscript_process$new(callr::rscript_process_options(
    script = files$script, cmdargs = c(files$job, files$result),
    stdout = files$out, stderr = files$err))
  list(process = process, files = files)
}

# This function returns the value of a finished job, or stops with its error
# and, when it died without one, the last lines it wrote.
.job_value <- function(job, label) {
  res <- if (file.exists(job$files$result)) readRDS(job$files$result)
  if (is.null(res)) {
    tail <- utils::tail(c(.read_lines(job$files$out), .read_lines(job$files$err)), 10L)
    .befa_stop(label, " stopped without a result. Its last messages:\n",
               paste(tail, collapse = "\n"))
  }
  if (!isTRUE(res$ok)) .befa_stop(label, " failed: ", res$message)
  res$value
}

# This function stops a job if it still runs and removes its files.
.stop_job <- function(job) {
  if (job$process$is_alive()) job$process$kill()
  unlink(unlist(job$files))
}

# This function returns the largest iteration a job has reported, from the
# lines "Iteration: N" that Stan writes while it samples.
.reported_iteration <- function(job) {
  lines <- c(.read_lines(job$files$out), .read_lines(job$files$err))
  at <- regmatches(lines, regexpr("Iteration:\\s*[0-9]+", lines))
  if (!length(at)) return(0L)
  max(as.integer(sub("\\D+", "", at)))
}

# This function reads a file that another process may still be writing.
.read_lines <- function(file) {
  if (!file.exists(file)) return(character(0))
  tryCatch(readLines(file, warn = FALSE), error = function(e) character(0))
}

# This function writes a duration for people.
.format_elapsed <- function(secs) {
  if (secs < 60) return(sprintf("%.1fs", secs))
  sprintf("%dm %ds", as.integer(secs %/% 60), as.integer(round(secs %% 60)))
}

# This function returns the style of the bar: light blue squares for the
# iterations run and a grey line for those to come. Their shapes differ, so
# the bar reads on light and dark consoles and without colours; a console that
# cannot draw them keeps cli's plain bar (NULL).
.bar_style <- function() {
  if (!cli::is_utf8_output()) return(NULL)
  done <- cli::make_ansi_style("#6CB4EE")("\u25a0")
  list(complete = done, current = done,
       incomplete = cli::make_ansi_style("#9AA0A6")("\u2500"))
}
