# Every test uses its own cache folder and a fixed backend status.
local_cache <- function(env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  withr::local_options(BayesianEFA.cache_dir = dir, .local_envir = env)
  local_mocked_bindings(.backend_status = function(...) {
    lapply(.backend_table(), function(x) list(available = TRUE, version = "1.0", reason = NULL))
  }, .env = env)
  dir
}

test_that("a program is named by the hash of its code", {
  a <- .code_hash("parameters { real y; }")
  expect_identical(a, .code_hash("parameters { real y; }"))
  expect_false(identical(a, .code_hash("parameters { real z; }")))
  expect_match(a, "^[0-9a-f]{32}$")
})

test_that("a program is written to the cache once, under its backend and version", {
  cache <- local_cache()
  code <- .stan_code(.befa_spec("cor", "unit_vector", 2L))
  files <- .model_files(code, "cmdstanr")
  expect_identical(normalizePath(dirname(files$dir)), normalizePath(cache))
  expect_match(basename(files$dir), "^befa-cmdstanr-1.0-[0-9a-f]{32}$")
  expect_identical(paste(readLines(files$stan_file), collapse = "\n"), sub("\n$", "", code))
  old <- as.POSIXct("2020-01-01 12:00:00")
  Sys.setFileTime(files$stan_file, old)
  .model_files(code, "cmdstanr")
  expect_lt(abs(as.numeric(file.mtime(files$stan_file)) - as.numeric(old)), 2)
})

test_that("each backend leaves its own compiled file", {
  local_cache()
  files <- .model_files("model {}", "rstan")
  expect_false(.is_compiled(files))
  file.create(file.path(files$dir, "model.rds"))
  expect_true(.is_compiled(files))

  files <- .model_files("model {}", "cmdstanr")
  file.create(file.path(files$dir, "model.exe"))
  expect_true(.is_compiled(files))
  files <- .model_files("model {}", "stanr")
  expect_false(.is_compiled(files))
  file.create(file.path(files$dir, ".model.so"))
  expect_true(.is_compiled(files))
  expect_true(.is_compiled(.model_files("model {}", "stanli")))
})

test_that("a model compiles only the first time it is needed", {
  local_cache()
  compiled <- 0L
  local_mocked_bindings(.run_jobs = function(jobs, ...) {
    compiled <<- compiled + 1L
    file.create(file.path(dirname(jobs[[1L]]$args$stan_file), "model.rds"))
    list(TRUE)
  })
  code <- .stan_code(.befa_spec("cov", "normal", 2L))
  expect_message(files <- .compiled_model(code, "rstan"), "once per model")
  .compiled_model(code, "rstan")
  expect_identical(compiled, 1L)
  .compiled_model(code, "stanli", quiet = TRUE)
  expect_identical(compiled, 1L)
})

test_that("only the most recently used models are kept", {
  cache <- local_cache()
  dirs <- file.path(cache, paste0("befa-x-", 1:4))
  for (i in seq_along(dirs)) {
    dir.create(dirs[i])
    Sys.setFileTime(dirs[i], as.POSIXct("2020-01-01") + i * 3600)
  }
  dir.create(file.path(cache, "something-else"))
  expect_length(.prune_cache(keep = 10L), 0L)
  removed <- .prune_cache(keep = 2L)
  expect_setequal(basename(removed), c("befa-x-1", "befa-x-2"))
  expect_setequal(basename(.cached_models()), c("befa-x-3", "befa-x-4"))
  expect_true(dir.exists(file.path(cache, "something-else")))
})

test_that("the cache can be emptied", {
  cache <- local_cache()
  expect_identical(befa_clear_cache(), 0L)
  dir.create(file.path(cache, "befa-x-1"))
  expect_identical(befa_clear_cache(), 1L)
  expect_length(.cached_models(), 0L)
  withr::local_options(BayesianEFA.cache_dir = file.path(cache, "nowhere"))
  expect_length(.cached_models(), 0L)
})

test_that("every structure can be compiled in advance", {
  local_cache()
  local_mocked_bindings(.resolve_backend = function(backend, ...) backend)
  seen <- character(0)
  local_mocked_bindings(.compiled_model = function(code, backend, quiet) {
    seen <<- c(seen, code)
    list(dir = "a-folder")
  })
  dirs <- befa_compile_models("cmdstanr", quiet = TRUE)
  expect_length(dirs, length(.spec_grid()))
  expect_identical(length(unique(seen)), 20L)
  expect_message(none <- befa_compile_models("stanli"), "needs no compilation")
  expect_length(none, 0L)
})

test_that("the grid holds the 16 models and the 4 null models", {
  grid <- .spec_grid()
  expect_length(grid, 20L)
  expect_identical(sum(vapply(grid, function(s) s$loadings == "none", TRUE)), 4L)
})
