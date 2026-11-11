# bridge-sampling.R
# Marginal likelihood by bridge sampling on a canonical representative of
# the rotations of the loadings
# Author: Ricardo Rey-Sáez
# Last modified: 11-11-2026

# Public functions

#' Marginal Likelihood by Bridge Sampling
#'
#' Estimates the log marginal likelihood (evidence) of a fitted model with
#' bridge sampling, taking into account that its loadings are not identified:
#' the posterior is invariant under every rotation, reflection and
#' permutation of the factors.
#'
#' @param object A `befa` object returned by [befa()], preferably fitted by
#'   MCMC, that did not draw from the prior alone (`prior_only = TRUE`).
#' @param method Bridge sampling method of [bridgesampling::bridge_sampler()]:
#'   `"warp3"` (default) or `"normal"`.
#' @param repetitions Number of repetitions of the estimate, which give its
#'   Monte Carlo error.
#' @param anchors Optional indices of the `n_factors` items whose rows define
#'   the canonical axes (see Details). By default they are chosen from the
#'   draws.
#' @param ... Other arguments of [bridgesampling::bridge_sampler()], such as
#'   `maxiter`.
#'
#' @details
#' Because the posterior is constant along the orbits of the orthogonal group
#' \eqn{O(M)}, MCMC chains visit only part of the equivalent copies of each
#' configuration, and bridge sampling applied to the raw draws is biased and
#' highly variable. The integral over the group is instead taken exactly:
#' with the QR decomposition of the block of `n_factors` anchor rows, every
#' draw of the loadings is written as \eqn{Z = TQ}, with \eqn{T} lower
#' trapezoidal with a positive diagonal and \eqn{Q \in O(M)}, so that
#' \deqn{p(y) = \mathrm{Vol}(O(M)) \int q(T)\,J(T)\,dT,\qquad
#' J(T) = \prod_{i=1}^{M} t_{ii}^{M-i},}
#' and bridge sampling estimates the integral over \eqn{T} from the
#' canonicalized draws, without sampling any new model. The anchors are
#' chosen by a pivoted Cholesky decomposition of the posterior mean of
#' \eqn{\Lambda\Lambda^\top} (strong items of distinct factors), with half of
#' the draws. With the unit-vector prior, the rows of \eqn{T} are written in
#' polar coordinates, since the prior on their norms places them on thin
#' spherical shells; with the normal prior and with one factor, in Cartesian
#' coordinates. With one factor, \eqn{O(1) = \{-1, 1\}} and the identity is
#' \eqn{p(y) = 2 p_+(y)}.
#'
#' The log density is computed in C++ with the formulas of the Stan
#' program, every constant included. The identity needs a prior on the
#' loadings that every rotation leaves unchanged: the unit-vector prior, the
#' normal prior with location 0 (its default) and, with one factor, a
#' symmetric Beta prior (equal shapes).
#'
#' @return The `bridge` object of [bridgesampling::bridge_sampler()] (a
#'   `bridge_list` with repetitions), whose `logml` already includes
#'   \eqn{\log \mathrm{Vol}(O(M))}, so it works with [bridgesampling::bf()]
#'   and [bridgesampling::post_prob()]. Its `mcse_logml` is the Monte Carlo
#'   standard error of `logml` of Micaletto and Vehtari (2025), reliable
#'   below 0.3. Its attribute `"befa"` holds the anchors, their pivots
#'   (residual communalities) and the coordinates used.
#'
#' @references
#' Gronau, Q. F., Singmann, H., & Wagenmakers, E.-J. (2020). bridgesampling:
#' An R package for estimating normalizing constants. *Journal of
#' Statistical Software, 92*(10), 1-29. \doi{10.18637/jss.v092.i10}
#'
#' Micaletto, G., & Vehtari, A. (2025). *Bridge sampling diagnostics*.
#' arXiv:2508.14487. \doi{10.48550/arXiv.2508.14487}
#'
#' Meng, X.-L., & Wong, W. H. (1996). Simulating ratios of normalizing
#' constants via a simple identity: A theoretical exploration. *Statistica
#' Sinica, 6*(4), 831-860.
#'
#' Muirhead, R. J. (1982). *Aspects of multivariate statistical theory*.
#' Wiley.
#'
#' @examples
#' \dontrun{
#' fit2 <- befa(HS_data, n_factors = 2, seed = 17, iter_sampling = 5000)
#' fit3 <- befa(HS_data, n_factors = 3, seed = 17, iter_sampling = 5000)
#' ml2 <- befa_bridge_sampler(fit2)
#' ml3 <- befa_bridge_sampler(fit3)
#' bridgesampling::bf(ml3, ml2)
#' }
#'
#' @export
befa_bridge_sampler <- function(object, method = c("warp3", "normal"),
                                repetitions = 1, anchors = NULL, ...) {
  .check_befa_object(object)
  if (.is_prior_only(object)) {
    .befa_stop("Bridge sampling needs draws from the posterior; this fit drew ",
               "from the prior alone (`prior_only = TRUE`).")
  }
  method <- match.arg(method)
  for (pkg in c("bridgesampling", "coda")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      .befa_stop("Bridge sampling needs the ", pkg, " package; install it with ",
                 "install.packages(\"", pkg, "\").")
    }
  }
  if (object$method != "mcmc") {
    .befa_warn("The draws come from ", .method_name(object$method), ", not ",
               "from MCMC, so the estimate of the marginal likelihood may be ",
               "biased.")
  }
  pars <- .stan_parameter_draws(object)
  setup <- .bridge_setup(object$spec, object$stan_data, pars,
                         n_chains = dim(object$draws)[2L],
                         G = .loading_gram(object), anchors = anchors)
  out <- bridgesampling::bridge_sampler(
    samples = setup$samples, log_posterior = .bridge_log_posterior,
    data = setup$data, lb = setup$lb, ub = setup$ub, method = method,
    repetitions = repetitions, cores = 1, silent = TRUE, ...)
  attr(out, "befa") <- setup$info
  out
}

# Internal functions

# This function returns the log density of the coordinates of bridge
# sampling (the function bridge_sampler() calls).
.bridge_log_posterior <- function(pars, data) {
  bridge_log_posterior_cpp(data$density, data$map, unname(pars))
}

# This function prepares bridge sampling from the draws of the Stan
# parameters (`pars`, one matrix per parameter, draws in rows, chains
# stacked): it canonicalizes the loadings of every draw, writes them in the
# coordinates of bridge sampling, checks on a few draws that the density is
# invariant and that the coordinates map back, and returns the samples (one
# coda::mcmc per chain), the data of the log density, the bounds and a
# summary. `G` is the posterior mean of Lambda Lambda', from which the
# anchors are chosen; `prior_only` drops the likelihood, so the marginal
# likelihood is 1 (for the tests).
.bridge_setup <- function(spec, stan_data, pars, n_chains, G, anchors = NULL,
                          prior_only = FALSE) {
  if (spec$loadings == "normal" && any(stan_data$pr_lambda_loc != 0)) {
    .befa_stop("Bridge sampling needs a prior on the loadings that rotations ",
               "leave unchanged; the normal prior must have location 0.")
  }
  if (spec$loadings == "uni" && stan_data$pr_lambda[1L] != stan_data$pr_lambda[2L]) {
    .befa_stop("Bridge sampling needs a prior on the loadings that a change of ",
               "sign leaves unchanged; the Beta prior must have equal shapes.")
  }
  J <- stan_data$J
  M <- if (spec$loadings == "uni") 1L else as.integer(stan_data$M)
  U <- .unconstrain(spec, pars)
  n_u <- ncol(U)
  load_idx <- switch(spec$loadings, uv = J + seq_len(J * M),
                     uni = seq_len(J), normal = seq_len(J * M))
  other_idx <- setdiff(seq_len(n_u), load_idx)
  polar <- spec$loadings == "uv"
  S <- nrow(U)
  first <- rep(seq_len(S / n_chains) <= round(S / n_chains / 2), n_chains)

  # Anchors and the canonical loadings of every draw
  piv <- .select_anchors(G, M, anchors)
  if (min(piv$pivots) < 0.05) {
    .befa_warn("Weak anchors (residual communality ", round(min(piv$pivots), 3),
               "): the estimate stays valid but may be less precise.")
  }
  ord <- c(piv$anchors, setdiff(seq_len(J), piv$anchors))
  T_list <- lapply(seq_len(S), function(s) {
    .rotate_to_T(matrix(U[s, load_idx], J, M), ord, M)
  })

  # Coordinates of T: polar rows need the chart centres of the rows past M
  centers <- bases <- rep(list(numeric(0)), J)
  if (polar) {
    for (i in seq_len(J - M) + M) {
      m <- Reduce(`+`, lapply(T_list[first], function(T) T[i, ] / sqrt(sum(T[i, ]^2))))
      centers[[i]] <- if (sum(m^2) > 1e-12) m / sqrt(sum(m^2)) else diag(M)[, 1L]
      bases[[i]] <- qr.Q(qr(cbind(centers[[i]], diag(M))))[, -1L, drop = FALSE]
    }
    phi <- t(vapply(T_list, .T_to_polar, numeric(J * M - M * (M - 1) / 2),
                    M = M, centers = centers, bases = bases))
  } else {
    phi <- t(vapply(T_list, .T_to_cartesian, numeric(J * M - M * (M - 1) / 2),
                    M = M))
  }

  log_vol <- .log_vol_orthogonal(M)
  density <- .befa_density(spec, stan_data, prior_only)
  map <- bridge_map_cpp(list(type = if (polar) 0L else 1L, p = J, M = M, n_u = n_u,
                             inv_ord = order(ord), load_idx = load_idx,
                             other_idx = other_idx, centers = centers,
                             bases = bases, log_vol = log_vol))

  # Checks on a few draws: the coordinates map back to T, and the density
  # does not change from Z to its canonical T
  for (s in seq_len(min(5L, S))) {
    back <- bridge_coords_to_T_cpp(map, phi[s, ])
    if (max(abs(back - T_list[[s]])) > 1e-8) {
      .befa_stop("The coordinates of bridge sampling do not map back to T.")
    }
    u_T <- U[s, ]
    u_T[load_idx] <- T_list[[s]][order(ord), , drop = FALSE]
    lp <- befa_log_density_cpp(density, U[s, ])
    if (abs(lp - befa_log_density_cpp(density, u_T)) > 1e-6 * max(1, abs(lp))) {
      .befa_stop("The log density is not invariant under rotations of the ",
                 "loadings, so the quotient does not apply.")
    }
  }

  samples <- cbind(U[, other_idx, drop = FALSE], phi)
  colnames(samples) <- c(sprintf("u[%d]", other_idx),
                         sprintf("phi[%d]", seq_len(ncol(phi))))
  chain <- rep(seq_len(n_chains), each = S / n_chains)
  samples <- coda::mcmc.list(lapply(seq_len(n_chains), function(ch) {
    coda::mcmc(samples[chain == ch, , drop = FALSE])
  }))
  names_all <- colnames(samples[[1L]])
  list(samples = samples,
       data = list(density = density, map = map),
       lb = stats::setNames(rep(-Inf, length(names_all)), names_all),
       ub = stats::setNames(rep(Inf, length(names_all)), names_all),
       info = list(coordinates = if (polar) "polar" else "cartesian", M = M,
                   anchors = piv$anchors, pivots = piv$pivots,
                   log_vol = log_vol))
}

# This function returns the draws of the parameters of the Stan program of a
# fit (draws in rows, chains stacked), which befa() keeps in `$fit$draws`.
.stan_parameter_draws <- function(object) {
  draws <- object$fit$draws
  names <- c(switch(object$spec$loadings, uv = c("h2", "Z"), uni = "lambda_std",
                    normal = c("Lambda", "psi")),
             if (.spec_has_scale(object$spec)) "sigma",
             if (object$model_type == "raw") "nu")
  stored <- unique(sub("\\[.*", "", posterior::variables(draws)))
  if (!all(names %in% stored)) {
    .befa_stop("The fit does not keep the draws of its Stan parameters; fit ",
               "the model again with this version of BayesianEFA.")
  }
  stats::setNames(lapply(names, function(n) .stacked_draws(draws, n)), names)
}

# This function returns the posterior mean of Lambda Lambda' (standardized),
# which the rotations leave unchanged, from the first half of each chain.
.loading_gram <- function(object) {
  L <- .stacked_draws(.metric_draws(object, "std"), "Lambda_std")
  dims <- dim(object$draws)
  first <- rep(seq_len(dims[1L]) <= round(dims[1L] / 2), dims[2L])
  J <- object$stan_data$J
  M <- object$n_factors
  G <- matrix(0, J, J)
  for (s in which(first)) G <- G + tcrossprod(matrix(L[s, ], J, M))
  G / sum(first)
}

# This function returns the unconstrained coordinates of draws of the Stan
# parameters (a list with one matrix per parameter, draws in rows), in the
# order of the parameters of the Stan program: logit of h2, log of sigma and
# psi, the logit of (lambda_std + 1) / 2, and Z, Lambda and nu as they are.
.unconstrain <- function(spec, pars) {
  parts <- switch(spec$loadings,
    uv = list(stats::qlogis(pars$h2), pars$Z),
    uni = list(stats::qlogis((pars$lambda_std + 1) / 2)),
    normal = list(pars$Lambda, log(pars$psi))
  )
  if (.spec_has_scale(spec)) parts <- c(parts, list(log(pars$sigma)))
  if (spec$model_type == "raw") parts <- c(parts, list(pars$nu))
  unname(do.call(cbind, parts))
}

# This function builds the C++ log density of a model structure from its
# Stan data: the priors, and the patterns of observed items in the form of
# the vectorized Woodbury likelihood.
.befa_density <- function(spec, stan_data, prior_only = FALSE) {
  sd <- stan_data
  J <- sd$J
  M <- if (spec$loadings == "uni") 1L else as.integer(sd$M)
  raw <- spec$model_type == "raw"
  pieces <- .stan_data_patterns(sd, raw)
  P <- length(pieces)
  obs_mat <- matrix(0, P, J)
  ybar_pad <- matrix(0, P, J)
  for (p in seq_len(P)) {
    obs_mat[p, pieces[[p]]$obs] <- 1
    if (raw) ybar_pad[p, pieces[[p]]$obs] <- pieces[[p]]$y_bar
  }
  root_pad <- do.call(rbind, lapply(pieces, function(pc) {
    out <- matrix(0, ncol(pc$root), J)
    out[, pc$obs] <- t(pc$root)
    out
  }))
  n_p <- vapply(pieces, `[[`, 0, "n")
  data <- list(obs_mat = obs_mat, root_pad = root_pad,
               n_obs = drop(crossprod(obs_mat, n_p)),
               r2_item = colSums(root_pad^2), ybar_pad = ybar_pad, n_p = n_p,
               K_p = vapply(pieces, function(pc) ncol(pc$root), 0),
               ll_const = -0.5 * sum(n_p * rowSums(obs_mat)) * log(2 * pi))

  get <- function(x, default = 0) if (is.null(x)) default else x
  xi <- get(sd$pr_xi)
  model <- list(
    loadings = match(spec$loadings, c("uv", "uni", "normal")) - 1L,
    scale = .spec_has_scale(spec), raw = raw, prior_only = prior_only,
    J = J, M = M,
    h2_a = get(sd$pr_h2[1L], 1), h2_b = get(sd$pr_h2[2L], 1), xi = xi,
    log_C = 0.5 * M * log(pi) + 0.5 * (xi + M) * log(2) +
      lgamma(0.5 * (xi + M)) - lgamma(0.5 * M),
    lam_a = get(sd$pr_lambda[1L], 1), lam_b = get(sd$pr_lambda[2L], 1),
    lam_loc = get(sd$pr_lambda_loc, rep(0, J)),
    lam_scale = get(sd$pr_lambda_scale, rep(1, J)),
    psi_shape = get(sd$pr_psi_shape, 1),
    psi_scale = get(sd$pr_psi_scale, rep(1, J)),
    sig_df = get(sd$pr_sigma[1L], 1), sig_loc = get(sd$pr_sigma[2L]),
    sig_scale = get(sd$pr_sigma[3L], 1),
    nu_loc = get(sd$pr_nu[1L]), nu_scale = get(sd$pr_nu[2L], 1)
  )
  befa_density_cpp(model, data)
}

# This function returns the patterns of the Stan data of a model: for each
# one, its rows, its observed items, its scatter root and (raw model) its
# means.
.stan_data_patterns <- function(sd, raw) {
  if (is.null(sd$P)) {
    return(list(list(n = sd$N, obs = seq_len(sd$J), root = sd$S_root,
                     y_bar = if (raw) sd$y_bar)))
  }
  pos_j <- 0L
  pos_r <- 0L
  lapply(seq_len(sd$P), function(p) {
    obs <- sd$obs_idx[pos_j + seq_len(sd$J_p[p])]
    root <- matrix(sd$root_flat[pos_r + seq_len(sd$J_p[p] * sd$K_p[p])],
                   sd$J_p[p], sd$K_p[p])
    y_bar <- if (raw) sd$ybar_flat[pos_j + seq_len(sd$J_p[p])]
    pos_j <<- pos_j + sd$J_p[p]
    pos_r <<- pos_r + sd$J_p[p] * sd$K_p[p]
    list(n = sd$n_p[p], obs = obs, root = root, y_bar = y_bar)
  })
}

# This function returns log Vol(O(M)), the volume of the orthogonal group
# (Muirhead, 1982, Theorem 2.1.15).
.log_vol_orthogonal <- function(M) {
  M * log(2) + M * (M + 1) / 4 * log(pi) - sum(lgamma(seq_len(M) / 2))
}

# This function chooses M anchor items with a greedy pivoted Cholesky
# decomposition of G = E[Lambda Lambda']: each anchor is the item with the
# largest communality left unexplained by the previous ones. It returns the
# anchors (or checks those given) and their pivots.
.select_anchors <- function(G, M, fixed = NULL) {
  J <- nrow(G)
  if (!is.null(fixed) &&
      (length(fixed) != M || anyDuplicated(fixed) || any(!fixed %in% seq_len(J)))) {
    .befa_stop("`anchors` must hold ", M, " distinct item indices.")
  }
  anchors <- integer(0)
  pivots <- numeric(0)
  for (m in seq_len(M)) {
    d <- diag(G)
    d[anchors] <- -Inf
    a <- if (is.null(fixed)) which.max(d) else fixed[m]
    anchors <- c(anchors, a)
    pivots <- c(pivots, G[a, a])
    G <- G - tcrossprod(G[, a]) / G[a, a]
  }
  list(anchors = as.integer(unname(anchors)), pivots = pivots)
}

# This function returns the canonical representative T of a draw of the
# loadings (Z[ord, ] = T Q, T lower trapezoidal with a positive diagonal).
.rotate_to_T <- function(Z, ord, M) {
  Zp <- Z[ord, , drop = FALSE]
  if (M == 1L) return(Zp * sign(Zp[1L, 1L]))
  qrA <- qr(t(Zp[seq_len(M), , drop = FALSE]))
  if (qrA$rank < M || any(qrA$pivot != seq_len(M))) {
    .befa_stop("The block of anchor rows of a draw is singular.")
  }
  sweep(Zp %*% qr.Q(qrA), 2L, sign(diag(qr.R(qrA))), `*`)
}

# This function writes T in Cartesian coordinates: the log of the diagonal,
# the entries below it in the anchor rows, and the rows past M.
.T_to_cartesian <- function(T, M) {
  lower <- unlist(lapply(seq_len(M)[-1L], function(i) T[i, seq_len(i - 1L)]))
  c(log(diag(T)[seq_len(M)]), lower, as.vector(t(T[-seq_len(M), , drop = FALSE])))
}

# This function writes T in polar coordinates: the log norm of every row,
# then a chart of the direction of the anchor rows 2..M (on hemispheres,
# around the pole) and of the rows past M (on full spheres, around their
# mean direction). A chart maps the angle theta from its centre to
# z = 2 atanh(theta / theta_max) v, almost linear for concentrated directions.
.T_to_polar <- function(T, M, centers, bases) {
  p <- nrow(T)
  r <- sqrt(rowSums(T^2))
  tau <- T / r
  z_hemi <- unlist(lapply(seq_len(M)[-1L], function(i) {
    w <- tau[i, seq_len(i - 1L)]
    nw <- sqrt(sum(w^2))
    theta <- atan2(nw, tau[i, i])
    if (nw > 0) 2 * atanh(2 * theta / pi) * w / nw else 0 * w
  }))
  z_sph <- unlist(lapply(seq_len(p - M) + M, function(i) {
    w <- drop(crossprod(bases[[i]], tau[i, ]))
    nw <- sqrt(sum(w^2))
    theta <- atan2(nw, sum(centers[[i]] * tau[i, ]))
    if (nw > 0) 2 * atanh(theta / pi) * w / nw else 0 * w
  }))
  c(log(r), z_hemi, z_sph)
}
