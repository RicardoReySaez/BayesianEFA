# backend-registry.R
# Which Stan backends the package can use and how one is chosen
# Author: Ricardo Rey-Sáez
# Last modified: 31-10-2026

# The backends never run in the R session of the user: they compile and sample
# in background processes. rstan, for instance, fails to load on Windows in a
# session where stanr is loaded, and a session that only orchestrates also
# stays free of Stan's messages.

# Internal functions

# The backends in the order `backend = "auto"` tries them, from the benchmark
# in bench/ (see bench/README.md). stanli needs no compilation, and on models
# up to 30 items it samples as fast as or faster than the compiled backends
# once those are cached; with 60 items, cmdstanr is 13% faster per step. The
# compiled backends follow, with rstan last for its older Stan.
.backend_order <- c("stanli", "cmdstanr", "stanr", "rstan")

# The estimation methods: Hamiltonian Monte Carlo (NUTS), automatic
# differentiation variational inference and Pathfinder.
.befa_methods <- c("mcmc", "advi", "pathfinder")

# This function returns what the package needs from each backend: its R
# package and minimum version, the methods it runs, and whether it compiles
# each model (stanli interprets it instead).
.backend_table <- function() {
  all <- .befa_methods
  list(
    rstan = list(package = "rstan", version = "2.26.0",
                 methods = c("mcmc", "advi"), compiles = TRUE),
    cmdstanr = list(package = "cmdstanr", version = "0.9.0",
                    methods = all, compiles = TRUE),
    stanr = list(package = "stanr", version = "0.1.5",
                 methods = all, compiles = TRUE),
    stanli = list(package = "stanli", version = "0.25.0",
                  methods = all, compiles = FALSE)
  )
}

# Backend checks found in this session, so each check runs once.
.befa_env <- new.env(parent = emptyenv())

# This function returns, for every backend, whether it can run and why not.
# The packages are found without loading them; what needs a loaded package
# (the CmdStan installation, the stanli runtime) is checked in a background
# session, once per R session.
.backend_status <- function(refresh = FALSE) {
  if (!refresh && !is.null(.befa_env$backend_status)) {
    return(.befa_env$backend_status)
  }
  table <- .backend_table()
  status <- lapply(names(table), function(b) {
    info <- table[[b]]
    version <- .installed_version(info$package)
    reason <- if (is.null(version)) {
      paste0("the ", info$package, " package is not installed")
    } else if (version < info$version) {
      paste0(info$package, " ", version, " is older than the ", info$version,
             " the package needs")
    } else if (b == "rstan") {
      .rstan_mismatch(version, .installed_version("StanHeaders"))
    }
    list(available = is.null(reason), version = as.character(version),
         reason = reason)
  })
  names(status) <- names(table)

  # Requirements that only a loaded backend can report
  candidates <- names(status)[vapply(status, `[[`, TRUE, "available")]
  extra <- .probe_backends(candidates)
  for (b in names(extra)) {
    if (!is.null(extra[[b]])) {
      status[[b]]$available <- FALSE
      status[[b]]$reason <- extra[[b]]
    }
  }
  .befa_env$backend_status <- status
  status
}

# This function tells why rstan cannot compile: a missing package it compiles
# against, or StanHeaders of another Stan release, since rstan only compiles
# against the headers of its own major and minor version (rstan 2.32 with
# StanHeaders 2.39 fails inside its C++ code). It returns NULL when rstan can
# compile.
.rstan_mismatch <- function(rstan, headers,
                            installed = function(p) nzchar(system.file(package = p))) {
  needs <- c("StanHeaders", "BH", "RcppEigen", "RcppParallel")
  absent <- needs[!vapply(needs, installed, TRUE)]
  if (length(absent)) {
    return(paste0("rstan compiles against ", paste(absent, collapse = ", "),
                  ", which ", if (length(absent) > 1L) "are" else "is",
                  " not installed"))
  }
  release <- function(v) paste(unlist(v)[1:2], collapse = ".")
  if (release(rstan) == release(headers)) return(NULL)
  paste0("rstan ", rstan, " cannot compile with StanHeaders ", headers,
         "; install versions of the same Stan release, for instance both ",
         "from https://stan-dev.r-universe.dev")
}

# This function returns the installed version of a package, or NULL, without
# loading it.
.installed_version <- function(package) {
  if (!nzchar(system.file(package = package))) return(NULL)
  utils::packageVersion(package)
}

# This function checks, in a background session, what an installed backend
# still needs: CmdStan for cmdstanr and the runtime for stanli. It returns,
# per backend, NULL or the reason it cannot run.
.probe_backends <- function(backends) {
  if (!length(backends)) return(list())
  callr::r(function(backends) {
    check <- function(b) {
      switch(b,
        cmdstanr = if (is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
          "cmdstanr cannot find CmdStan; install it with cmdstanr::install_cmdstan()"
        },
        stanli = if (!isTRUE(stanli::stanli_available())) {
          "the stanli runtime is not installed; install it with stanli::stanli_install()"
        },
        NULL)
    }
    out <- lapply(backends, function(b) tryCatch(check(b), error = function(e) {
      paste0(b, " failed to load: ", conditionMessage(e))
    }))
    names(out) <- backends
    out
  }, args = list(backends = backends))
}

# This function chooses the backend: the one requested, if it can run the
# method, or with `backend = "auto"` the first one in .backend_order that is
# installed and runs the method.
.resolve_backend <- function(backend = "auto", method = "mcmc",
                             verbose = TRUE) {
  method <- .check_method(method)
  table <- .backend_table()
  if (!identical(backend, "auto")) {
    if (!is.character(backend) || length(backend) != 1L ||
        !backend %in% names(table)) {
      .befa_stop("`backend` must be \"auto\" or one of ",
                 paste0("\"", names(table), "\"", collapse = ", "), ".")
    }
    if (!method %in% table[[backend]]$methods) {
      .befa_stop("The ", backend, " backend does not run ", .method_name(method),
                 "; use another backend, such as \"cmdstanr\" or \"stanli\".")
    }
    status <- .backend_status()[[backend]]
    if (!status$available) {
      .befa_stop("The ", backend, " backend cannot run: ", status$reason, ".")
    }
    return(backend)
  }

  status <- .backend_status()
  usable <- .backend_order[vapply(.backend_order, function(b) {
    status[[b]]$available && method %in% table[[b]]$methods
  }, TRUE)]
  if (!length(usable)) {
    .befa_stop("No Stan backend can run ", .method_name(method), ". Install ",
               "one of them, for instance stanli, which needs no compiler:\n",
               "  install.packages(\"stanli\", repos = ",
               "\"https://seantalts.r-universe.dev\"); stanli::stanli_install()")
  }
  if (verbose) .befa_msg("using the ", usable[1L], " backend.")
  usable[1L]
}

# This function checks the estimation method.
.check_method <- function(method) {
  if (!is.character(method) || length(method) != 1L ||
      !method %in% .befa_methods) {
    .befa_stop("`method` must be one of ",
               paste0("\"", .befa_methods, "\"", collapse = ", "), ".")
  }
  method
}

# This function names a method in messages.
.method_name <- function(method) {
  c(mcmc = "MCMC", advi = "ADVI", pathfinder = "Pathfinder")[[method]]
}
