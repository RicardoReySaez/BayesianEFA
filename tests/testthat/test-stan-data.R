test_that("rows are grouped by the items they observe, complete rows first", {
  Y <- rbind(c(1, NA, 3), c(1, 2, 3), c(NA, 2, NA), c(4, NA, 6), c(1, 2, 3))
  groups <- .missing_patterns(Y)
  expect_identical(groups[[1L]], list(rows = c(2L, 5L), obs = 1:3))
  expect_identical(lapply(groups[-1L], `[[`, "rows"), list(3L, c(1L, 4L)))
  expect_identical(groups[[3L]]$obs, c(1L, 3L))
})

test_that("listwise deletion keeps the complete rows and says how many it removed", {
  Y <- .simulated_data(missing = TRUE)
  expect_warning(prep <- .prepare_befa_data(Y, "raw", "listwise"), "rows with missing")
  expect_identical(prep$N, sum(stats::complete.cases(Y)))
  expect_identical(.data_patterns(prep), "single")
})

test_that("FIML drops only the rows without any observed value", {
  Y <- .simulated_data(missing = TRUE)
  Y[3L, ] <- NA
  expect_warning(prep <- .prepare_befa_data(Y, "cov", "FIML"), "without any observed")
  expect_identical(prep$N, nrow(Y) - 1L)
  expect_identical(.data_patterns(prep), "multiple")
  expect_error(suppressWarnings(.prepare_befa_data(matrix(NA_real_, 3, 2), "cov", "FIML")),
               "No rows are left")
})

test_that("the correlation model standardizes and the covariance model centers", {
  Y <- .simulated_data()
  cor_prep <- .prepare_befa_data(Y, "cor")
  # Standardizing with the ML SDs makes the cross-products N times cor(Y).
  expect_equal(crossprod(cor_prep$Y), nrow(Y) * stats::cor(Y), tolerance = 1e-12)
  cov_prep <- .prepare_befa_data(Y, "cov")
  expect_equal(colMeans(cov_prep$Y), rep(0, ncol(Y)), tolerance = 1e-12)
  expect_equal(.prepare_befa_data(Y, "raw")$Y, unname(Y))
})

test_that("summary statistics give one pattern with the ML matrix", {
  Y <- .simulated_data()
  S <- stats::cov(Y)
  prep <- .prepare_befa_data(model_type = "raw", sample_nobs = 40, sample_cov = S,
                             sample_mean = colMeans(Y))
  expect_true(prep$summary)
  expect_equal(prep$S, S * 39 / 40)
  expect_identical(.data_patterns(prep), "single")
  cor_prep <- .prepare_befa_data(model_type = "cor", sample_nobs = 40, sample_cov = S)
  expect_equal(cor_prep$S, unname(stats::cov2cor(S)))
  cor_prep2 <- .prepare_befa_data(model_type = "cor", sample_nobs = 40,
                                  sample_cor = stats::cor(Y))
  expect_equal(cor_prep2$S, cor_prep$S)
})

test_that("summary statistics are checked", {
  S <- diag(3)
  expect_error(.prepare_befa_data(model_type = "cor", sample_cor = S), "sample_nobs")
  expect_error(.prepare_befa_data(model_type = "cov", sample_nobs = 10, sample_cor = S),
               "needs `sample_cov`")
  expect_error(.prepare_befa_data(model_type = "raw", sample_nobs = 10, sample_cov = S),
               "needs `sample_mean`")
  expect_error(.prepare_befa_data(model_type = "cov", sample_nobs = 10,
                                  sample_cov = "a"), "numeric matrix")
})

test_that("scatter roots rebuild the scatter matrix with few columns", {
  set.seed(51)
  many <- matrix(stats::rnorm(30), 10, 3)
  root <- .scatter_root(many)
  expect_identical(dim(root), c(3L, 3L))
  expect_equal(tcrossprod(root), crossprod(many), tolerance = 1e-12)

  few <- matrix(stats::rnorm(6), 2, 3)
  expect_identical(.scatter_root(few), t(few))

  # Two identical columns make the scatter matrix singular.
  x <- stats::rnorm(8)
  singular <- unname(cbind(x, x, stats::rnorm(8)))
  root <- .scatter_root(singular)
  expect_identical(ncol(root), 2L)
  expect_equal(tcrossprod(root), crossprod(singular), tolerance = 1e-10)

  expect_identical(.matrix_root(matrix(0, 2, 2)), matrix(0, 2, 1))

  # A Cholesky factor that ends with a pivot of rounding size, as some LAPACK
  # libraries return for a singular matrix, is not used. Adding 1e-14 along
  # the null direction lets chol() succeed on every platform.
  S <- crossprod(singular)
  v <- c(1, -1, 0) / sqrt(2)
  nearly <- S + 1e-14 * max(diag(S)) * tcrossprod(v)
  expect_false(is.null(tryCatch(chol(nearly), error = function(e) NULL)))
  root <- .matrix_root(nearly)
  expect_identical(ncol(root), 2L)
  expect_equal(tcrossprod(root), nearly, tolerance = 1e-10)
})

test_that("item SDs come from the data, the EM or the summary statistics", {
  Y <- .simulated_data()
  expect_equal(.item_sd(.prepare_befa_data(Y, "raw")), apply(Y, 2, .ml_sd))
  Y_na <- .simulated_data(missing = TRUE)
  prep <- .prepare_befa_data(Y_na, "cov", "FIML")
  expect_equal(.item_sd(prep), sqrt(diag(estimate_h1_moments_em(prep$Y)$Sigma)))
  sp <- .prepare_befa_data(model_type = "cov", sample_nobs = 40, sample_cov = stats::cov(Y))
  expect_equal(.item_sd(sp), sqrt(diag(sp$S)))
})

test_that("a single pattern passes one root and the means", {
  spec <- .befa_spec("raw", "unit_vector", 2L)
  prep <- .prepared_for(spec)
  sd <- .stan_data(spec, prep, 2L, .prior_for(spec))
  Yc <- sweep(prep$Y, 2, colMeans(prep$Y))
  expect_equal(tcrossprod(sd$S_root), crossprod(Yc), tolerance = 1e-10)
  expect_equal(sd$y_bar, colMeans(prep$Y))
  expect_identical(sd$K, ncol(sd$S_root))
  expect_identical(sd$M, 2L)
  expect_named(sd, c("N", "J", "M", "K", "S_root", "y_bar", "prior_only",
                     "pr_h2", "pr_xi", "pr_sigma", "pr_nu"))
})

test_that("the data switch the likelihood off only when asked", {
  spec <- .befa_spec("cor", "unit_vector", 2L)
  prep <- .prepared_for(spec)
  expect_identical(.stan_data(spec, prep, 2L, .prior_for(spec))$prior_only, 0L)
  expect_identical(.stan_data(spec, prep, 2L, .prior_for(spec),
                              prior_only = TRUE)$prior_only, 1L)
})

test_that("several patterns stack their items, roots and means", {
  spec <- .befa_spec("raw", "normal", 2L, patterns = "multiple")
  prep <- .prepared_for(spec)
  sd <- .stan_data(spec, prep, 2L, .prior_for(spec))
  expect_identical(sd$P, length(prep$groups))
  expect_identical(sum(sd$n_p), prep$N)
  expect_identical(sd$n_root, as.integer(sum(sd$J_p * sd$K_p)))
  expect_length(sd$obs_idx, sum(sd$J_p))
  expect_length(sd$ybar_flat, sum(sd$J_p))

  # Each stacked root rebuilds the scatter of its rows around their means.
  pos_j <- 0
  pos_r <- 0
  for (p in seq_len(sd$P)) {
    g <- prep$groups[[p]]
    Yg <- prep$Y[g$rows, g$obs, drop = FALSE]
    Yc <- sweep(Yg, 2, colMeans(Yg))
    R <- matrix(sd$root_flat[pos_r + seq_len(sd$J_p[p] * sd$K_p[p])], sd$J_p[p])
    expect_equal(tcrossprod(R), crossprod(Yc), tolerance = 1e-10)
    expect_equal(sd$obs_idx[pos_j + seq_len(sd$J_p[p])], g$obs)
    pos_j <- pos_j + sd$J_p[p]
    pos_r <- pos_r + sd$J_p[p] * sd$K_p[p]
  }
})

test_that("the structure and the data must agree on the patterns", {
  spec <- .befa_spec("cov", "unit_vector", 2L, patterns = "multiple")
  prep <- .prepare_befa_data(.simulated_data(), "cov")
  expect_error(.stan_data(spec, prep, 2L, .prior_for(spec)), "expects multiple")
})

test_that("the normal prior is scaled by the item SDs", {
  spec <- .befa_spec("cov", "normal", 2L)
  prior <- list(lambda = c(0.5, 2), psi = c(3, 4))
  item_sd <- c(1, 2, 3)
  pd <- .stan_prior_data(spec, prior, item_sd)
  expect_equal(pd$pr_lambda_loc, 0.5 * item_sd)
  expect_equal(pd$pr_lambda_scale, 2 * item_sd)
  expect_identical(pd$pr_psi_shape, 3)
  expect_equal(pd$pr_psi_scale, 4 * item_sd^2)
  expect_named(.stan_prior_data(.befa_spec("cor", "unit_vector", 1L), list(lambda = c(2, 2))),
               "pr_lambda")
})

test_that("the prepared data keep the names of the items", {
  Y <- .simulated_data()
  expect_identical(.prepare_befa_data(Y, "cor")$item_names, paste0("Item_", 1:5))
  colnames(Y) <- paste0("q", 1:5)
  expect_identical(.prepare_befa_data(Y, "cor")$item_names, paste0("q", 1:5))
  expect_identical(.prepare_befa_data(as.data.frame(Y), "raw", "FIML")$item_names,
                   paste0("q", 1:5))
  summ <- .prepare_befa_data(model_type = "raw", sample_nobs = 40,
                             sample_cov = unname(stats::cov(Y)),
                             sample_mean = colMeans(Y))
  expect_identical(summ$item_names, paste0("q", 1:5))
  expect_identical(.prepare_befa_data(model_type = "cor", sample_nobs = 40,
                                      sample_cor = stats::cor(Y))$item_names,
                   paste0("q", 1:5))
  # Incomplete or repeated names are replaced as a whole.
  expect_identical(.item_labels(c("a", "", "c"), 3), paste0("Item_", 1:3))
  expect_identical(.item_labels(c("a", "a", "c"), 3), paste0("Item_", 1:3))
  expect_identical(.item_labels(c("a", NA, "c"), 3), paste0("Item_", 1:3))
})
