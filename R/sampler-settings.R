# sampler-settings.R
# The options of each estimation method, with their defaults and checks
# Author: Ricardo Rey-Sáez
# Last modified: 31-10-2026

# Internal functions

# This function completes and checks the options of a method given through
# `...` of befa(). The names follow cmdstanr; for MCMC the rstan names `iter`
# (warmup included), `warmup`, `cores` and `control` are also accepted.
# Unknown options are errors, so a misspelled one is never ignored.
.sampler_settings <- function(method, dots = list()) {
  method <- .check_method(method)
  dots <- .rstan_aliases(method, dots)
  defaults <- .sampler_defaults(method)
  unknown <- setdiff(names(dots), names(defaults))
  if (length(unknown)) {
    .befa_stop("`", unknown[1L], "` is not an option of ", .method_name(method),
               ". It accepts ", paste0("`", names(defaults), "`", collapse = ", "),
               ".")
  }
  s <- utils::modifyList(defaults, dots)

  if (is.null(s$seed)) s$seed <- sample.int(.Machine$integer.max, 1L)
  .check_count(s$seed, "seed", lower = 0)
  s$seed <- as.integer(s$seed)

  if (method == "mcmc") {
    for (arg in c("chains", "parallel_chains", "iter_warmup", "iter_sampling",
                  "max_treedepth")) {
      .check_count(s[[arg]], arg)
      s[[arg]] <- as.integer(s[[arg]])
    }
    .check_probability(s$adapt_delta, "adapt_delta")
    total <- s$iter_warmup + s$iter_sampling
    if (is.null(s$refresh)) s$refresh <- max(1L, total %/% 100L)
  } else {
    for (arg in c("iter", "draws", if (method == "pathfinder") "num_paths")) {
      .check_count(s[[arg]], arg)
      s[[arg]] <- as.integer(s[[arg]])
    }
    if (method == "advi") {
      s$algorithm <- .check_choice(s$algorithm, c("meanfield", "fullrank"),
                                   "algorithm")
    }
    if (is.null(s$refresh)) s$refresh <- max(1L, s$iter %/% 100L)
  }
  if (!is.numeric(s$refresh) || length(s$refresh) != 1L || s$refresh < 0) {
    .befa_stop("`refresh` must be zero (no progress bar) or a positive number.")
  }
  s$refresh <- as.integer(s$refresh)
  s$method <- method
  s
}

# This function returns the default options of a method. MCMC runs four
# chains of 1000 warmup and 1000 kept iterations with max_treedepth = 15, as
# in v0.1.0, and adapt_delta = 0.9: in bench/bench-nuts.R, 0.95 gave the
# smallest effective sample size per second on every data set, and no value
# from 0.8 up gave a divergent transition. ADVI and Pathfinder return 1000
# draws; Pathfinder runs four paths and resamples them with PSIS.
.sampler_defaults <- function(method) {
  common <- list(seed = NULL, refresh = NULL, init = NULL)
  switch(method,
    mcmc = c(list(chains = 4L, parallel_chains = .default_cores(),
                  iter_warmup = 1000L, iter_sampling = 1000L,
                  adapt_delta = 0.9, max_treedepth = 15L), common),
    advi = c(list(algorithm = "meanfield", iter = 10000L, draws = 1000L,
                  tol_rel_obj = 0.01), common),
    pathfinder = c(list(iter = 1000L, draws = 1000L, num_paths = 4L,
                        psis_resample = TRUE), common)
  )
}

# This function translates the rstan names of the MCMC options: `iter`
# counts the warmup too, `cores` sets the chains that run at once and
# `control` holds adapt_delta and max_treedepth.
.rstan_aliases <- function(method, dots) {
  if (method != "mcmc") return(dots)
  if (!is.null(dots$iter)) {
    warmup <- if (is.null(dots$warmup)) dots$iter %/% 2 else dots$warmup
    if (dots$iter <= warmup) {
      .befa_stop("`iter` counts the warmup too, so it must exceed `warmup`.")
    }
    dots$iter_warmup <- warmup
    dots$iter_sampling <- dots$iter - warmup
    dots$iter <- NULL
    dots$warmup <- NULL
  } else if (!is.null(dots$warmup)) {
    dots$iter_warmup <- dots$warmup
    dots$warmup <- NULL
  }
  if (!is.null(dots$cores)) {
    dots$parallel_chains <- dots$cores
    dots$cores <- NULL
  }
  if (!is.null(dots$control)) {
    unknown <- setdiff(names(dots$control), c("adapt_delta", "max_treedepth"))
    if (length(unknown)) {
      .befa_stop("`control` accepts `adapt_delta` and `max_treedepth`, not `",
                 unknown[1L], "`.")
    }
    dots[names(dots$control)] <- dots$control
    dots$control <- NULL
  }
  dots
}

# This function returns how many chains run at once by default: up to four,
# leaving two cores free, and no more than two while R CMD check limits the
# cores a package may use.
.default_cores <- function() {
  n <- tryCatch(parallel::detectCores(), error = function(e) NA_integer_)
  cores <- if (is.na(n)) 1L else max(1L, min(n - 2L, 4L))
  limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(limit) && !identical(tolower(limit), "false")) {
    cores <- min(cores, 2L)
  }
  as.integer(cores)
}
