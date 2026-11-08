# befa-fit-measures.R
# Bayesian fit indices (DevM framework) and information criteria
# Author: Ricardo Rey-Sáez
# Last modified: 08-11-2026

# Public functions

#' Bayesian Fit Indices for Exploratory Factor Analysis
#'
#' Computes the Bayesian analogues of the usual structural equation modeling
#' fit indices proposed by Garnier-Villarreal and Jorgensen (2020), at every
#' posterior draw, together with information criteria of out-of-sample
#' predictive accuracy (PSIS-LOO, WAIC, DIC, DIC_i and BPIC).
#'
#' @param object A `befa` object returned by [befa()], fitted to raw `data`.
#' @param loo_args List with `r_eff` (logical; relative efficiencies of the
#'   MCMC draws) and `cores`, passed to [loo::loo()]. Defaults to the
#'   `loo_args` of the fit.
#'
#' @details
#' **_The DevM framework_**
#'
#' The deviance of the model at each draw, relative to the saturated model,
#' plays the role of the \eqn{\chi^2} statistic, \eqn{D_i}. The effective
#' number of parameters \eqn{pD} is the \eqn{p_{loo}} of PSIS-LOO (Vehtari et
#' al., 2017), the number of sample moments is \eqn{p^*} (\eqn{J(J+3)/2} with
#' means, \eqn{J(J+1)/2} for a covariance matrix and \eqn{J(J-1)/2} for a
#' correlation matrix), the Bayesian degrees of freedom are
#' \eqn{df_B = p^* - pD}, and the noncentrality is
#' \eqn{\lambda_i = \max(0, D_i - p^*)}.
#'
#' **_Indices_**
#'
#' * BRMSEA \eqn{= \sqrt{\lambda_i / (df_B N)}}; BGamma \eqn{= J / (J + 2
#'   \lambda_i / N)} and its adjusted version \eqn{1 - (p^* / df_B)(1 -
#'   \hat\Gamma)}; BMc \eqn{= \exp(-\lambda_i / (2N))}; and the SRMR of each
#'   draw.
#' * Chi2_ppp: the posterior predictive p-value, the share of draws whose
#'   replicated data deviate more than the observed data (ideally above .05).
#' * BCFI and BTLI compare the model with the null (independence) model, which
#'   is fitted with the same backend and method (in the covariance and raw
#'   models; the null correlation model has no free parameters). Both are
#'   bounded to \[0, 1\].
#'
#' **_Information criteria_**
#'
#' They use the marginal log-likelihood of each person (the factor scores
#' integrated out), with deviance \eqn{D = -2 \log p(y \mid \theta)}.
#'
#' * ELPD, LOOIC and p_loo from PSIS-LOO (Vehtari et al., 2017), and WAIC
#'   with \eqn{p_{WAIC}}, the sum of the posterior variances of the pointwise
#'   log-likelihood (Watanabe, 2010).
#' * DIC \eqn{= \bar D + p_D}, with \eqn{p_D = \bar D - D(\bar\theta)}
#'   (Spiegelhalter et al., 2002); \eqn{\bar\theta} holds the posterior means
#'   of the aligned parameters. A negative \eqn{p_D} means those means
#'   describe the posterior poorly (for instance, with `rotate = "none"`).
#' * DIC_i \eqn{= \bar D + p_V}, with \eqn{p_V = \mathrm{Var}(D) / 2}, the
#'   parameterization-invariant DIC of Xiao and Rabe-Hesketh (2026), which
#'   needs no point estimate.
#' * BPIC \eqn{= \bar D + 2 p_D} (Ando, 2007).
#'
#' None of them is a function of a single draw, so they have no posterior
#' distribution. Their standard errors are computed across persons: from the
#' pointwise contributions, or by the jackknife over persons for DIC_i and
#' \eqn{p_V}, whose variance term does not split into persons. With ADVI and
#' Pathfinder, PSIS-LOO treats the draws as independent, and the
#' approximation can bias every criterion.
#'
#' With missing values, the log-likelihood uses every observed value (FIML)
#' and the saturated model is estimated with the EM algorithm.
#'
#' @return A list of class `befa_fitmeasures` with:
#' * `$fit_draws`: A `draws_array` of `Chi2`, `Chi2_Null`, `BRMSEA`,
#'   `BGamma`, `Adj_BGamma`, `BMc`, `SRMR`, `BCFI` and `BTLI`.
#' * `$summary`: Their posterior summaries.
#' * `$posterior_fit`: The same draws as a data frame.
#' * `$information`: A matrix with the estimate and standard error of ELPD,
#'   LOOIC, p_loo, WAIC, p_WAIC, DIC, p_D, DIC_i, p_V and BPIC.
#' * `$loo_object`: The `loo` object of the model.
#' * `$details`: `p_star`, `pD`, `N` and `chi2_ppp`.
#'
#' @references
#' Ando, T. (2007). Bayesian predictive information criterion for the
#' evaluation of hierarchical Bayesian and empirical Bayes models.
#' *Biometrika, 94*(2), 443-458.
#'
#' Garnier-Villarreal, M., & Jorgensen, T. D. (2020). Adapting fit indices for
#' Bayesian structural equation modeling: Comparison to maximum likelihood.
#' *Psychological Methods, 25*(1), 46-70. \doi{10.1037/met0000224}
#'
#' Spiegelhalter, D. J., Best, N. G., Carlin, B. P., & van der Linde, A.
#' (2002). Bayesian measures of model complexity and fit. *Journal of the
#' Royal Statistical Society: Series B, 64*(4), 583-639.
#'
#' Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC. *Statistics and
#' Computing, 27*(5), 1413-1432. \doi{10.1007/s11222-016-9696-4}
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and
#' widely applicable information criterion in singular learning theory.
#' *Journal of Machine Learning Research, 11*, 3571-3594.
#'
#' Xiao, X., & Rabe-Hesketh, S. (2026). *A parameterization-invariant DIC*.
#' arXiv:2605.27844. \doi{10.48550/arXiv.2605.27844}
#'
#' @examples
#' \dontrun{
#' fit <- befa(HS_data, n_factors = 3, seed = 17)
#' befa_fit_measures(fit)
#' }
#'
#' @export
befa_fit_measures <- function(object, loo_args = NULL) {
  if (inherits(object, "befa_list")) {
    return(lapply(object, befa_fit_measures, loo_args = loo_args))
  }
  .check_befa_object(object)
  if (!is.null(object$fit_indices)) return(object$fit_indices)
  if (object$data$summary) {
    .befa_warn("Fit measures need the raw `data`; this model was fitted to ",
               "summary statistics, so they are not computed.")
    return(NULL)
  }
  .fit_measures(object, validate_loo_args(
    if (is.null(loo_args)) object$options$loo_args else loo_args))
}

#' Print Bayesian Fit Measures
#'
#' Displays the fit indices and the information criteria of a
#' `befa_fitmeasures` object as two tables.
#'
#' @param x A `befa_fitmeasures` object from [befa_fit_measures()].
#' @param digits Number of decimal places.
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#' @method print befa_fitmeasures
#' @export
print.befa_fitmeasures <- function(x, digits = 3, ...) {
  .print_fit_tables(build_fit_summary_table(x), build_information_table(x),
                    c(0.025, 0.975), digits, first = 1L)
  invisible(x)
}

# Internal functions

# This function returns a value held in the environment `shared`, computing
# and keeping it the first time. Fits with several numbers of factors share
# in this way what does not depend on the number of factors; without
# `shared`, the value is simply computed.
.shared_value <- function(shared, name, value) {
  if (is.null(shared)) return(value)
  if (is.null(shared[[name]])) shared[[name]] <- value
  shared[[name]]
}

# This function computes the fit indices and information criteria of a fit
# to raw data. `shared` holds the saturated and null models when several fits
# use the same data, so each of them is estimated only once.
.fit_measures <- function(object, loo_config, shared = NULL) {
  saturated <- .shared_value(shared, "saturated",
                             .saturated_moments(object$data$Y, object$model_type))
  metrics <- .posterior_metrics(object, saturated)
  loo_res <- .loo_fit(metrics$log_lik, object$method, loo_config)
  p_loo <- loo_res$estimates["p_loo", "Estimate"]
  information <- .information_criteria(object, metrics$log_lik, loo_res)
  null <- .shared_value(shared, "null",
                        .null_metrics(object, saturated$ll, loo_config))

  N <- nrow(object$data$Y)
  J <- object$stan_data$J
  p_star <- .sample_moments(object$model_type, J)

  # Absolute indices
  lambda_hat <- pmax(0, metrics$chisq - p_star)
  df_H <- max(1, p_star - p_loo)
  bgamma <- J / (J + (2 / N) * lambda_hat)

  # Incremental indices against the null model
  ratio_null <- null$lambda_vec / null$df
  ratio_H <- lambda_hat / df_H
  posterior_df <- data.frame(
    Chi2 = metrics$chisq,
    Chi2_Null = null$chisq_vec,
    BRMSEA = sqrt(lambda_hat / (df_H * N)),
    BGamma = bgamma,
    Adj_BGamma = 1 - (p_star / df_H) * (1 - bgamma),
    BMc = exp(-0.5 * lambda_hat / N),
    SRMR = metrics$srmr,
    BCFI = pmin(1, pmax(0, 1 - lambda_hat / null$lambda_vec)),
    BTLI = pmin(1, pmax(0, (ratio_null - ratio_H) / (ratio_null - 1)))
  )

  dims <- dim(object$draws)
  fit_draws <- .as_chains(as.matrix(posterior_df), dims[1L], dims[2L])
  structure(list(
    fit_draws = fit_draws,
    summary = posterior::summarise_draws(fit_draws),
    posterior_fit = posterior_df,
    information = information,
    loo_object = loo_res,
    details = list(p_star = p_star, pD = p_loo, N = N,
                   chi2_ppp = mean(metrics$chisq_rep >= metrics$chisq))
  ), class = "befa_fitmeasures")
}

# This function returns the number of sample moments the data model fits.
.sample_moments <- function(model_type, J) {
  switch(model_type, raw = J * (J + 3) / 2, cov = J * (J + 1) / 2,
         cor = J * (J - 1) / 2)
}

# This function returns the PSIS-LOO of a pointwise log-likelihood array
# (iterations x chains x persons). MCMC draws can use their relative
# efficiencies; the draws of ADVI and Pathfinder are independent.
.loo_fit <- function(log_lik, method, loo_config) {
  r_eff <- if (method == "mcmc" && isTRUE(loo_config$r_eff)) {
    loo::relative_eff(exp(log_lik), cores = loo_config$cores)
  } else if (method == "mcmc") {
    NA
  } else {
    1
  }
  suppressWarnings(loo::loo(log_lik, r_eff = r_eff, cores = loo_config$cores))
}

# This function returns the information criteria of a fit (rows) with their
# estimate and standard error across persons (columns), from the pointwise
# log-likelihood array (iterations x chains x persons): PSIS-LOO and WAIC
# from the loo package; DIC, p_D and BPIC summed person by person, with the
# deviance at the posterior means of the aligned parameters; and DIC_i and
# p_V, whose variance term does not split into persons, with the jackknife
# over persons. Leaving person i out adds 2 log p(y_i | theta) to the
# deviance of every draw, so the jackknife needs no refit.
.information_criteria <- function(object, log_lik, loo_res) {
  ll <- matrix(log_lik, prod(dim(log_lik)[1:2]))
  N <- ncol(ll)
  mean_ll <- colMeans(ll)
  p_d <- 2 * (.plugin_loglik(object) - mean_ll)
  waic <- suppressWarnings(loo::waic(ll))$estimates
  person_sum <- function(x) c(sum(x), sqrt(N * stats::var(x)))

  D <- -2 * rowSums(ll)
  p_v <- stats::var(D) / 2
  jack_p_v <- (stats::var(D) + 4 * apply(ll, 2L, stats::var) +
                 4 * drop(stats::cov(D, ll))) / 2
  jack_dic_i <- mean(D) + 2 * mean_ll + jack_p_v
  jack_se <- function(x) sqrt((N - 1) / N * sum((x - mean(x))^2))

  out <- rbind(
    ELPD = loo_res$estimates["elpd_loo", ],
    LOOIC = loo_res$estimates["looic", ],
    p_loo = loo_res$estimates["p_loo", ],
    WAIC = waic["waic", ],
    p_WAIC = waic["p_waic", ],
    DIC = person_sum(-2 * mean_ll + p_d),
    p_D = person_sum(p_d),
    DIC_i = c(mean(D) + p_v, jack_se(jack_dic_i)),
    p_V = c(p_v, jack_se(jack_p_v)),
    BPIC = person_sum(-2 * mean_ll + 2 * p_d)
  )
  colnames(out) <- c("Estimate", "SE")
  out
}

# This function returns the log-likelihood of each person at the posterior
# means of the aligned parameters, in the metric the model was fitted to,
# pattern by pattern of observed items.
.plugin_loglik <- function(object) {
  p <- .fitted_parameters(object)
  L <- matrix(colMeans(p$lambda), p$J, p$M)
  Sigma <- tcrossprod(L) + diag(colMeans(p$psi), p$J)
  nu <- colMeans(p$nu)
  Y <- object$data$Y
  out <- numeric(nrow(Y))
  for (g in object$data$groups) {
    o <- g$obs
    out[g$rows] <- mvnfast::dmvn(Y[g$rows, o, drop = FALSE], nu[o],
                                 Sigma[o, o, drop = FALSE], log = TRUE)
  }
  out
}

# This function returns the moments of the saturated model and its
# log-likelihood: the EM moments with missing values, the ML ones otherwise.
# They depend on the data and the data model, not on the number of factors.
.saturated_moments <- function(Y, model_type) {
  N <- nrow(Y)
  J <- ncol(Y)
  comp_idx <- which(stats::complete.cases(Y))
  miss_idx <- which(!stats::complete.cases(Y))
  if (length(miss_idx)) {
    h1 <- estimate_h1_moments_em(Y)
    mean <- if (model_type == "raw") h1$Mu else rep(0, J)
    cov <- if (model_type == "cor") stats::cov2cor(h1$Sigma) else h1$Sigma
    fiml <- list(comp_idx = comp_idx, miss_idx = miss_idx,
                 Y_comp = Y[comp_idx, , drop = FALSE],
                 Y_miss = Y[miss_idx, , drop = FALSE])
    ll <- sum(compute_fiml_loglik(fiml, mu = mean, Sigma = cov))
  } else {
    mean <- if (model_type == "raw") colMeans(Y) else rep(0, J)
    cov <- if (model_type == "cor") stats::cor(Y) else stats::cov(Y) * (N - 1) / N
    ll <- sum(mvnfast::dmvn(Y, mu = mean, sigma = cov, log = TRUE))
  }
  list(mean = mean, cov = cov, ll = ll)
}

# This function returns, for every draw, the pointwise log-likelihood of the
# rows the model was fitted to (FIML with missing values), the chi-square
# deviance against the saturated model, its replicated counterpart for the
# posterior predictive p-value, and the SRMR. The data are those of the fit,
# already standardized (correlation model) or centered (covariance model).
.posterior_metrics <- function(object,
                               saturated = .saturated_moments(object$data$Y,
                                                              object$model_type)) {
  Y <- object$data$Y
  N <- nrow(Y)
  J <- ncol(Y)
  comp_idx <- which(stats::complete.cases(Y))
  miss_idx <- which(!stats::complete.cases(Y))
  ll_saturated <- saturated$ll

  p <- .fitted_parameters(object)
  out <- compute_posterior_metrics_cpp(
    lambda = p$lambda, psi = p$psi, nu = p$nu, M = p$M, n_iter = p$n_iter,
    data = Y, M_sat = saturated$mean, S_sat = saturated$cov,
    ll_saturated = ll_saturated,
    model_type = match(object$model_type, c("raw", "cov", "cor")),
    has_missing = length(miss_idx) > 0L,
    Y_comp = if (length(comp_idx)) Y[comp_idx, , drop = FALSE] else matrix(0, 0, J),
    Y_miss = if (length(miss_idx)) Y[miss_idx, , drop = FALSE] else matrix(0, 0, J),
    comp_idx_0 = as.integer(comp_idx - 1L), miss_idx_0 = as.integer(miss_idx - 1L),
    obs_patterns_miss = lapply(miss_idx, function(i) as.integer(which(!is.na(Y[i, ])) - 1L))
  )
  c(out, list(ll_saturated = ll_saturated, N = N, J = J))
}

# This function returns the chi-square deviance of the null (independence)
# model at each draw, its noncentrality and its Bayesian degrees of freedom.
# The null covariance and raw models are fitted with the backend, method and
# settings of the fit (their programs compile once, like any other); the null
# correlation model has no free parameters, so its deviance is a constant.
.null_metrics <- function(object, ll_saturated, loo_config) {
  Y <- object$data$Y
  p_star <- .sample_moments(object$model_type, ncol(Y))
  dims <- dim(object$draws)

  if (object$model_type == "cor") {
    ll_null <- sum(stats::dnorm(Y[!is.na(Y)], log = TRUE))
    chisq <- 2 * (ll_saturated - ll_null)
    S <- dims[1L] * dims[2L]
    return(list(chisq_vec = rep(chisq, S), pD = 0, df = max(1, p_star),
                lambda_vec = rep(max(0, chisq - p_star), S)))
  }

  spec <- .befa_null_spec(object$model_type, object$spec$patterns)
  settings <- object$fit$settings
  settings$refresh <- 0L
  n <- switch(settings$method, mcmc = settings$chains, advi = 1L,
              pathfinder = settings$num_paths)
  inits <- .start_values(spec, object$data, object$n_factors,
                         object$priors_used, n, settings$seed)
  if (settings$method == "advi") inits <- inits[[1L]]
  fit <- .run_stan(.stan_code(spec), .stan_data(spec, object$data, NA, object$priors_used),
                   inits, object$backend, settings, quiet = TRUE)

  sigma <- .stacked_draws(fit$draws, "sigma")
  nu <- if (object$model_type == "raw") {
    .stacked_draws(fit$draws, "nu")
  } else {
    matrix(0, nrow(sigma), ncol(sigma))
  }
  ll <- null_loglik_cpp(Y, nu, sigma)
  chisq <- 2 * (ll_saturated - rowSums(ll))

  n_iter <- dim(fit$draws)[1L]
  log_lik <- array(NA_real_, dim = c(n_iter, dim(fit$draws)[2L], nrow(Y)))
  for (ch in seq_len(dim(fit$draws)[2L])) {
    log_lik[, ch, ] <- ll[(ch - 1L) * n_iter + seq_len(n_iter), , drop = FALSE]
  }
  pD <- .loo_fit(log_lik, settings$method, loo_config)$estimates["p_loo", "Estimate"]
  list(chisq_vec = chisq, pD = pD, df = max(1, p_star - pD),
       lambda_vec = pmax(0, chisq - p_star))
}
