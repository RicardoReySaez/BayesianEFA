# This function returns fake Stan draws (n_iter x n_chains) of the given
# parameters, named as Stan names them (column-major for matrices).
.fake_stan_draws <- function(pars, n_iter = 5L, n_chains = 2L) {
  cols <- lapply(names(pars), function(p) {
    x <- pars[[p]]
    if (is.matrix(x)) {
      names <- paste0(p, "[", rep(seq_len(nrow(x)), ncol(x)), ",",
                      rep(seq_len(ncol(x)), each = nrow(x)), "]")
      values <- c(x)
    } else {
      names <- paste0(p, "[", seq_along(x), "]")
      values <- x
    }
    m <- matrix(rep(values, each = n_iter * n_chains), nrow = n_iter * n_chains)
    colnames(m) <- names
    m
  })
  .as_chains(do.call(cbind, cols), n_iter, n_chains)
}

test_that("unit-vector loadings are rebuilt from Z and h2", {
  spec <- .befa_spec("cov", "unit_vector", 2L)
  Z <- matrix(c(3, 0, 4, 1, 2, 2), 3, 2)
  h2 <- c(0.5, 0.2, 0.8)
  draws <- .fake_stan_draws(list(h2 = h2, Z = Z, sigma = c(1, 2, 3)))
  std <- .standardized_parameters(draws, spec, 3L, 2L)
  L <- matrix(std$lambda[1, ], 3, 2)
  expect_equal(rowSums(L^2), h2)
  expect_equal(L[1, ] / sqrt(sum(L[1, ]^2)), Z[1, ] / sqrt(sum(Z[1, ]^2)))
  expect_equal(std$psi[1, ], 1 - h2)
  expect_equal(std$scale[1, ], c(1, 2, 3))
})

test_that("normal loadings are standardized by the implied item SDs", {
  spec <- .befa_spec("raw", "normal", 2L)
  Lambda <- matrix(c(1, 2, 0.5, 0, 1, 1), 3, 2)
  psi <- c(1, 2, 0.75)
  draws <- .fake_stan_draws(list(Lambda = Lambda, psi = psi, nu = c(0, 1, 2)))
  std <- .standardized_parameters(draws, spec, 3L, 2L)
  variance <- rowSums(Lambda^2) + psi
  expect_equal(std$scale[1, ], sqrt(variance))
  expect_equal(matrix(std$lambda[1, ], 3, 2), Lambda / sqrt(variance))
  expect_equal(std$psi[1, ], psi / variance)
  expect_equal(std$nu[1, ], c(0, 1, 2))
})

test_that("one factor under the unit-vector prior has uniquenesses 1 - lambda^2", {
  spec <- .befa_spec("cor", "unit_vector", 1L)
  draws <- .fake_stan_draws(list(lambda_std = c(0.5, -0.3)))
  std <- .standardized_parameters(draws, spec, 2L, 1L)
  expect_equal(std$psi[1, ], c(0.75, 0.91))
  expect_null(std$scale)
})

test_that("each model returns its parameters under the names of its metric", {
  final <- function(spec, M = 2L) {
    std <- list(lambda = matrix(0.5, 4, 3 * M), psi = matrix(0.4, 4, 3),
                scale = matrix(2, 4, 3), nu = matrix(1, 4, 3))
    if (spec$model_type != "raw") std$nu <- NULL
    colnames(.model_metric(std$lambda, std, spec, 3L, M))
  }
  cor <- final(.befa_spec("cor", "unit_vector", 2L))
  expect_identical(unique(sub("\\[.*", "", cor)), c("Lambda_std", "Psi_std", "h2"))
  raw_uv <- final(.befa_spec("raw", "unit_vector", 2L))
  expect_identical(unique(sub("\\[.*", "", raw_uv)), c("Lambda", "Psi", "sigma", "h2", "nu"))
  cov_normal <- final(.befa_spec("cov", "normal", 2L))
  expect_identical(unique(sub("\\[.*", "", cov_normal)), c("Lambda", "Psi"))
  expect_true("Lambda[3,2]" %in% cov_normal)
  uni <- final(.befa_spec("cor", "unit_vector", 1L), M = 1L)
  expect_true(all(c("Lambda_std[3,1]", "Psi_std[3]") %in% uni))
})

test_that("the final draws keep their chains and the metric scales the loadings", {
  spec <- .befa_spec("cov", "unit_vector", 2L)
  set.seed(61)
  Z <- matrix(stats::rnorm(6, 0, 3), 3, 2)
  draws <- .fake_stan_draws(list(h2 = c(0.5, 0.6, 0.4), Z = Z, sigma = c(1, 2, 3)))
  fit <- list(draws = draws)
  post <- .befa_posterior(fit, spec, 3L, 2L, rotate = "none")
  expect_null(post$rsp)
  # Lambda (6), Psi (3), sigma (3) and h2 (3)
  expect_identical(unname(dim(post$draws)), c(5L, 2L, 15L))
  L <- matrix(.stacked_draws(post$draws, "Lambda")[1, ], 3, 2)
  std <- matrix(.standardized_parameters(draws, spec, 3L, 2L)$lambda[1, ], 3, 2)
  expect_equal(L, std * c(1, 2, 3))
  expect_equal(as.numeric(.stacked_draws(post$draws, "Psi")[1, ]),
               c(0.5, 0.4, 0.6) * c(1, 4, 9))
})

test_that("the draws are aligned unless the rotation is switched off", {
  skip_on_cran()
  spec <- .befa_spec("cor", "unit_vector", 2L)
  set.seed(62)
  n <- 40L
  Z <- matrix(c(4, 4, 4, 0.5, 0.5, 0.5, 0.5, 0.5, 4, 4), 5, 2)
  # Each draw flips the sign of its first column, which the alignment undoes.
  rows <- lapply(seq_len(n), function(s) {
    z <- Z + matrix(stats::rnorm(10, 0, 0.1), 5, 2)
    if (s %% 2 == 0) z[, 1] <- -z[, 1]
    c(rep(0.6, 5), c(z))
  })
  m <- do.call(rbind, rows)
  colnames(m) <- c(paste0("h2[", 1:5, "]"),
                   paste0("Z[", rep(1:5, 2), ",", rep(1:2, each = 5), "]"))
  fit <- list(draws = .as_chains(m, n / 2, 2L))
  # The fake draws barely vary, which posterior notes while summarising them.
  post <- suppressWarnings(.befa_posterior(fit, spec, 5L, 2L))
  L <- .stacked_draws(post$draws, "Lambda_std")
  expect_named(post$rsp, c("Lambda_star", "objective", "sign_vectors", "perm_vectors"))
  expect_lt(max(apply(L, 2, stats::sd)), 0.1)
})

test_that("draws are converted to the standardized metric on demand", {
  spec <- .befa_spec("raw", "normal", 2L)
  Lambda <- matrix(c(1, 2, 0.5, 0, 1, 1), 3, 2)
  psi <- c(1, 2, 0.75)
  fit <- list(draws = .fake_stan_draws(list(Lambda = Lambda, psi = psi, nu = 1:3)))
  post <- .befa_posterior(fit, spec, 3L, 2L, rotate = "none")
  object <- list(draws = post$draws, model_type = "raw", n_factors = 2L,
                 stan_data = list(J = 3L))
  std <- .metric_draws(object, "std")
  variance <- rowSums(Lambda^2) + psi
  expect_equal(matrix(.stacked_draws(std, "Lambda_std")[1, ], 3, 2), Lambda / sqrt(variance))
  expect_equal(as.numeric(.stacked_draws(std, "h2")[1, ]), 1 - psi / variance)
  expect_identical(.metric_draws(object, "raw"), post$draws)
  expect_identical(.metric_draws(object), post$draws)

  object$model_type <- "cor"
  expect_identical(.metric_draws(object, "std"), post$draws)
  expect_error(.metric_draws(object, "raw"), "no metric of the items")
})
