# This function returns a backend status where only `available` can run.
.fake_status <- function(available) {
  out <- lapply(names(.backend_table()), function(b) {
    list(available = b %in% available, version = "9.9.9",
         reason = if (!b %in% available) paste(b, "is missing"))
  })
  names(out) <- names(.backend_table())
  out
}

test_that("every backend runs MCMC and ADVI, and rstan has no Pathfinder", {
  table <- .backend_table()
  expect_setequal(names(table), .backend_order)
  for (b in names(table)) expect_true(all(c("mcmc", "advi") %in% table[[b]]$methods))
  expect_false("pathfinder" %in% table$rstan$methods)
  expect_false(table$stanli$compiles)
})

test_that("rstan must compile against StanHeaders of its own Stan release", {
  all_there <- function(p) TRUE
  expect_null(.rstan_mismatch(package_version("2.32.7"), package_version("2.32.2"),
                              all_there))
  expect_match(.rstan_mismatch(package_version("2.32.7"), package_version("2.39.1"),
                               all_there), "cannot compile with StanHeaders 2.39.1")
  no_bh <- function(p) p != "BH"
  expect_match(.rstan_mismatch(package_version("2.32.7"), NULL, no_bh),
               "compiles against BH, which is not installed")
})

test_that("installed versions are read without loading the package", {
  expect_s3_class(.installed_version("stats"), "package_version")
  expect_null(.installed_version("no.such.package.here"))
})

test_that("the status of the backends is computed once per session", {
  calls <- 0L
  local_mocked_bindings(
    .installed_version = function(package) package_version("100.0.0"),
    .probe_backends = function(backends) {
      calls <<- calls + 1L
      list(stanli = "the stanli runtime is not installed")
    }
  )
  withr::defer(rm("backend_status", envir = .befa_env))
  status <- .backend_status(refresh = TRUE)
  expect_true(status$cmdstanr$available)
  expect_false(status$stanli$available)
  expect_match(status$stanli$reason, "runtime")
  .backend_status()
  expect_identical(calls, 1L)
})

test_that("old or missing backends are reported as such", {
  local_mocked_bindings(
    .installed_version = function(package) {
      if (package == "stanr") NULL else package_version("0.0.1")
    },
    .probe_backends = function(backends) list()
  )
  withr::defer(rm("backend_status", envir = .befa_env))
  status <- .backend_status(refresh = TRUE)
  expect_match(status$stanr$reason, "not installed")
  expect_match(status$cmdstanr$reason, "older than")
})

test_that("auto chooses the first usable backend for the method", {
  local_mocked_bindings(.backend_status = function(...) .fake_status(c("stanli", "rstan")))
  expect_message(b <- .resolve_backend("auto", "mcmc"), "using the stanli backend")
  expect_identical(b, "stanli")
  expect_identical(.resolve_backend("auto", "advi", verbose = FALSE), "stanli")
  local_mocked_bindings(.backend_status = function(...) .fake_status(c("cmdstanr", "stanli")))
  expect_identical(.resolve_backend("auto", "mcmc", verbose = FALSE), "stanli")

  local_mocked_bindings(.backend_status = function(...) .fake_status("rstan"))
  expect_identical(.resolve_backend("auto", "mcmc", verbose = FALSE), "rstan")
  expect_error(.resolve_backend("auto", "pathfinder"), "No Stan backend can run Pathfinder")
})

test_that("a requested backend must exist, run the method and be installed", {
  local_mocked_bindings(.backend_status = function(...) .fake_status("stanli"))
  expect_identical(.resolve_backend("stanli", "pathfinder"), "stanli")
  expect_error(.resolve_backend("brms"), "must be \"auto\" or one of")
  expect_error(.resolve_backend(c("stanli", "rstan")), "must be \"auto\"")
  expect_error(.resolve_backend("rstan", "pathfinder"), "does not run Pathfinder")
  expect_error(.resolve_backend("cmdstanr"), "cannot run: cmdstanr is missing")
})

test_that("methods are checked and named", {
  expect_error(.check_method("laplace"), "`method` must be one of")
  expect_identical(.method_name("advi"), "ADVI")
})

test_that("the background probe reports what an installed backend lacks", {
  skip_on_cran()
  expect_identical(.probe_backends(character(0)), list())
  skip_if_not_installed("stanli", "0.25.0")
  res <- .probe_backends("stanli")
  expect_named(res, "stanli")
  expect_true(is.null(res$stanli) || is.character(res$stanli))
})
