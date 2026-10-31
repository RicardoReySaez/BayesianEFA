# stan-cache.R
# Keep each compiled Stan model, so every model compiles once
# Author: Ricardo Rey-Sáez
# Last modified: 31-10-2026

# Public functions

#' Compile the Stan models of BayesianEFA in advance
#'
#' Every structure of the model (correlation, covariance or raw model, prior
#' on the loadings, complete or incomplete data) has its own Stan program,
#' which [befa()] compiles the first time it is needed and keeps in a cache
#' for later sessions. This function compiles all of them at once, for
#' instance right after installing the package. The stanli backend interprets
#' the programs and never compiles them.
#'
#' @param backend Stan backend: `"auto"` (the first one installed), or one of
#'   `"cmdstanr"`, `"stanr"` and `"rstan"`.
#' @param quiet Logical. If `TRUE`, no messages are shown.
#'
#' @return Invisibly, the folders of the compiled models.
#' @seealso [befa_clear_cache()] to remove them.
#' @export
befa_compile_models <- function(backend = "auto", quiet = FALSE) {
  backend <- .resolve_backend(backend, "mcmc", verbose = !quiet)
  if (!.backend_table()[[backend]]$compiles) {
    if (!quiet) .befa_msg("the ", backend, " backend needs no compilation.")
    return(invisible(character(0)))
  }
  specs <- .spec_grid()
  dirs <- vapply(specs, function(spec) {
    .compiled_model(.stan_code(spec), backend, quiet = quiet)$dir
  }, "")
  invisible(unname(dirs))
}

#' Remove the compiled Stan models of BayesianEFA
#'
#' Deletes the cache where [befa()] keeps the compiled Stan models. They are
#' compiled again the next time they are needed.
#'
#' @return Invisibly, the number of models removed.
#' @export
befa_clear_cache <- function() {
  models <- .cached_models()
  unlink(models, recursive = TRUE)
  invisible(length(models))
}

# Internal functions

# This function returns the folder of the cache. The option
# `BayesianEFA.cache_dir` moves it, for instance to a temporary folder in the
# tests.
.cache_dir <- function() {
  getOption("BayesianEFA.cache_dir", tools::R_user_dir("BayesianEFA", "cache"))
}

# This function returns the folders of the cached models.
.cached_models <- function() {
  cache <- .cache_dir()
  if (!dir.exists(cache)) return(character(0))
  entries <- list.files(cache, full.names = TRUE)
  entries[dir.exists(entries) & startsWith(basename(entries), "befa-")]
}

# This function returns the MD5 hash of a Stan program, which names its
# folder in the cache.
.code_hash <- function(code) {
  file <- tempfile(fileext = ".stan")
  on.exit(unlink(file), add = TRUE)
  writeLines(code, file, sep = "")
  unname(tools::md5sum(file))
}

# This function returns where a program lives in the cache for one backend:
# a folder named after the backend, its version and the hash of the code,
# which holds model.stan and what the backend compiles from it. The program
# is written only once, so its date stays that of its first compilation.
.model_files <- function(code, backend) {
  version <- .backend_status()[[backend]]$version
  key <- paste("befa", backend, version, .code_hash(code), sep = "-")
  dir <- file.path(.cache_dir(), key)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  stan_file <- file.path(dir, "model.stan")
  if (!file.exists(stan_file)) writeLines(code, stan_file, sep = "")
  list(key = key, dir = dir, stan_file = stan_file, backend = backend)
}

# This function tells whether the cache already holds what a backend compiles
# from a program: the rstan model object, the CmdStan executable or the
# shared library of stanr.
.is_compiled <- function(files) {
  listed <- list.files(files$dir, all.files = TRUE, no.. = TRUE)
  switch(files$backend,
    rstan = "model.rds" %in% listed,
    cmdstanr = any(listed %in% c("model", "model.exe")),
    stanr = any(grepl("^\\.model\\.(so|dll)$", listed)),
    stanli = TRUE
  )
}

# This function returns the cached model of a program, compiling it once if
# needed. Compilation runs in a background session, like sampling, behind a
# spinner. Using a model refreshes the date of its folder, and only the most
# recently used models are kept.
.compiled_model <- function(code, backend, quiet = FALSE) {
  files <- .model_files(code, backend)
  if (!.is_compiled(files)) {
    if (!quiet) {
      cli::cli_alert_info(paste("Compiling the Stan model with {backend};",
                                "this happens once per model."))
    }
    job <- list(func = .compile_job,
                args = list(backend = backend, stan_file = files$stan_file))
    .run_jobs(list(job), label = "Compiling the Stan model",
              done = "Compiled the Stan model", quiet = quiet)
  }
  Sys.setFileTime(files$dir, Sys.time())
  .prune_cache()
  files
}

# This function compiles a program with one backend, in the background
# session that runs it. It must not call anything of BayesianEFA, which that
# session may not have loaded.
.compile_job <- function(backend, stan_file) {
  switch(backend,
    rstan = {
      rstan::rstan_options(auto_write = TRUE)
      invisible(rstan::stan_model(file = stan_file, auto_write = TRUE))
    },
    cmdstanr = invisible(cmdstanr::cmdstan_model(stan_file, quiet = TRUE)),
    stanr = invisible(stanr::stan_model(stan_file = stan_file, quiet = TRUE))
  )
  TRUE
}

# This function keeps the `keep` most recently used models in the cache and
# removes the others. There are about 30 programs for each backend.
.prune_cache <- function(keep = getOption("BayesianEFA.cache_keep", 100L)) {
  models <- .cached_models()
  if (length(models) <= keep) return(invisible(character(0)))
  old <- models[order(file.mtime(models), decreasing = TRUE)][-seq_len(keep)]
  unlink(old, recursive = TRUE)
  invisible(old)
}

# This function lists every structure the generator supports: the 16 models
# and the 4 null models.
.spec_grid <- function() {
  specs <- list()
  for (model_type in c("cor", "cov", "raw")) {
    for (loadings in c("uv", "uni", "normal")) {
      if (model_type == "cor" && loadings == "normal") next
      for (patterns in c("single", "multiple")) {
        spec <- .befa_spec(
          model_type,
          lambda_prior = if (loadings == "normal") "normal" else "unit_vector",
          n_factors = if (loadings == "uni") 1L else 2L,
          patterns = patterns
        )
        specs[[.spec_label(spec)]] <- spec
      }
    }
  }
  for (model_type in c("cov", "raw")) {
    for (patterns in c("single", "multiple")) {
      spec <- .befa_null_spec(model_type, patterns)
      specs[[.spec_label(spec)]] <- spec
    }
  }
  specs
}
