# befa-compare.R
# Compare fitted models by fit, LOO, information criteria and model evidence
# Author: Ricardo Rey-Sáez
# Last modified: 08-11-2026

# Public functions

#' Compare BEFA Models
#'
#' Compares two or more models fitted to the same data, for instance with
#' different numbers of factors, in four tables, each ordered from the best
#' model to the worst: the Bayesian fit indices, the PSIS-LOO comparison with
#' stacking and pseudo-BMA+ weights, the information criteria with Akaike
#' weights, and the model evidence with posterior model probabilities.
#'
#' @param ... Two or more `befa` objects fitted to the same data with the
#'   same data model, or one list of them. Their names label the models
#'   (`Model 1`, `Model 2`, ... when they have none). A `befa_list` given
#'   alone returns the comparison it holds, if any; with other arguments, or
#'   without one, the comparison is computed.
#' @param criteria Tables to compute: any of `"fit"`, `"loo"`, `"ic"` and
#'   `"evidence"` (default all four). Fits to summary statistics have no
#'   pointwise log-likelihood, so they allow only `"evidence"`.
#' @param model_prior Prior model probabilities, one per model, for the
#'   posterior model probabilities. `NULL` (default) makes the models equally
#'   probable a priori.
#' @param repetitions Repetitions of bridge sampling per model (default 1).
#'   The Monte Carlo error of each evidence comes from the closed form of
#'   Micaletto and Vehtari (2025), so repetitions are not needed for it.
#' @param pseudobma_bb Logical. Use the Bayesian bootstrap for the pseudo-BMA
#'   weights (pseudo-BMA+, default) or not.
#' @param seed Seed of the Bayesian bootstrap and of the simulation of the
#'   posterior model probabilities, for reproducible results.
#' @param bridge_args List of other arguments of [befa_bridge_sampler()],
#'   such as `method` or `maxiter`.
#'
#' @details
#' **_Bayesian fit indices_**
#'
#' The posterior means of the fit indices of [befa_fit_measures()] and the
#' posterior predictive p-value, ordered by BRMSEA. A fit without
#' `$fit_indices` gets them computed, which fits its null model.
#'
#' **_LOO comparison_**
#'
#' Each model's ELPD from PSIS-LOO, ordered from the largest, its difference
#' from the best model (`elpd_diff`, with its standard error `se_diff`), and
#' `p_worse`, the probability that the model predicts worse than the best
#' one, \eqn{\Phi(-\mathrm{elpd\_diff} / \mathrm{se\_diff})}. The weights come
#' from [loo::loo_model_weights()]: stacking chooses the mixture of the
#' models' predictive distributions that predicts best (Yao et al., 2018),
#' and pseudo-BMA+ normalizes \eqn{\exp(\mathrm{elpd}_k)} over a Bayesian
#' bootstrap of the observations. Neither is a probability of a model.
#'
#' **_Information criteria_**
#'
#' WAIC, DIC, DIC_i and BPIC of each model (see [befa_fit_measures()]),
#' ordered by WAIC; lower values indicate better expected predictive
#' accuracy. Each comes with its Akaike weights,
#' \eqn{w_k = \exp(-\Delta_k / 2) / \sum_j \exp(-\Delta_j / 2)}, where
#' \eqn{\Delta_k} is the difference from the smallest value. LOOIC is
#' \eqn{-2} times the ELPD of the LOO comparison.
#'
#' **_Model evidence_**
#'
#' The log marginal likelihood of each model from [befa_bridge_sampler()],
#' ordered from the largest, with the Monte Carlo standard error of
#' Micaletto and Vehtari (2025) (with repetitions, their median), which is
#' reliable below 0.3; and the posterior model probabilities
#' \deqn{P(M_k \mid y) = \frac{p(y \mid M_k) P(M_k)}
#' {\sum_j p(y \mid M_j) P(M_j)}.}
#' Their intervals come from simulating the log evidences from their
#' estimates and Monte Carlo errors. Unlike the LOO weights, these
#' probabilities depend on the priors of every parameter, and assume that
#' one of the compared models generated the data.
#'
#' The fit indices, the LOO and the information criteria come from
#' `$fit_indices` when the fit holds them; otherwise they are computed.
#'
#' @return An object of class `befa_comparison` with `$models` (the names),
#'   `$fit` (the table of fit indices), `$loo` (the `loo` objects and the
#'   comparison table), `$ic` (the criteria and their Akaike weights),
#'   `$evidence` (the `bridge` objects, the table and the prior model
#'   probabilities) and `$notes`. Every table is ordered from the best model.
#'
#' @references
#' Micaletto, G., & Vehtari, A. (2025). *Bridge sampling diagnostics*.
#' arXiv:2508.14487. \doi{10.48550/arXiv.2508.14487}
#'
#' Yao, Y., Vehtari, A., Simpson, D., & Gelman, A. (2018). Using stacking to
#' average Bayesian predictive distributions. *Bayesian Analysis, 13*(3),
#' 917-1007. \doi{10.1214/17-BA1091}
#'
#' @examples
#' \dontrun{
#' fit2 <- befa(HS_data, n_factors = 2, seed = 17)
#' fit3 <- befa(HS_data, n_factors = 3, seed = 17)
#' befa_compare(two = fit2, three = fit3, seed = 1)
#' }
#'
#' @export
befa_compare <- function(..., criteria = c("fit", "loo", "ic", "evidence"),
                         model_prior = NULL, repetitions = 1L,
                         pseudobma_bb = TRUE, seed = NULL, bridge_args = list()) {
  models <- list(...)
  if (length(models) == 1L && is.list(models[[1L]]) && !inherits(models[[1L]], "befa")) {
    stored <- attr(models[[1L]], "comparison")
    if (nargs() == 1L && !is.null(stored)) return(stored)
    models <- models[[1L]]
  }
  if (length(models) < 2L || !all(vapply(models, inherits, TRUE, "befa"))) {
    .befa_stop("Compare two or more fits returned by befa().")
  }
  names(models) <- .model_names(names(models), length(models))
  criteria <- unique(match.arg(criteria, several.ok = TRUE))
  .check_count(repetitions, "repetitions")
  .check_flag(pseudobma_bb, "pseudobma_bb")
  .check_same_data(models)
  prior <- .model_prior(model_prior, length(models))

  notes <- character(0)
  if (models[[1L]]$data$summary && any(c("fit", "loo", "ic") %in% criteria)) {
    notes <- c(notes, paste("Fits to summary statistics have no pointwise",
                            "log-likelihood, so only the evidence is compared."))
    criteria <- intersect(criteria, "evidence")
  }

  out <- list(models = names(models), criteria = criteria)
  measures <- if (any(c("fit", "loo", "ic") %in% criteria)) {
    lapply(models, .compare_measures, fit = "fit" %in% criteria)
  }
  run <- function(expr) if (is.null(seed)) expr else withr::with_seed(seed, expr)
  if ("fit" %in% criteria) out$fit <- .compare_fit(measures)
  if ("loo" %in% criteria) {
    out$loo <- run(.compare_loo(lapply(measures, `[[`, "loo"), pseudobma_bb))
  }
  if ("ic" %in% criteria) out$ic <- .compare_ic(measures)
  if ("evidence" %in% criteria) {
    evidence <- .compare_evidence(models, repetitions, bridge_args)
    notes <- c(notes, evidence$notes)
    if (all(is.finite(evidence$log_ml))) {
      probs <- run(.posterior_model_probs(evidence$log_ml, evidence$mc_se, prior))
      tab <- cbind(log_ML = evidence$log_ml, MC_SE = evidence$mc_se, probs)
      rownames(tab) <- names(models)
      evidence$table <- tab[order(-tab[, "log_ML"]), , drop = FALSE]
    }
    evidence$model_prior <- stats::setNames(prior, names(models))
    out$evidence <- evidence[c("bridges", "table", "model_prior")]
  }
  out$pseudobma_bb <- pseudobma_bb
  out$notes <- notes
  structure(out, class = "befa_comparison")
}

#' Print a Comparison of BEFA Models
#'
#' Prints the tables of [befa_compare()], each ordered from the best model.
#'
#' @param x A `befa_comparison` object.
#' @param digits Number of decimal places of the ELPD, the criteria and the
#'   evidence (fit indices, weights and probabilities have three).
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#' @method print befa_comparison
#' @export
print.befa_comparison <- function(x, digits = 2, ...) {
  .print_comparison(x, digits)
  invisible(x)
}

# Internal functions

# This function prints the tables of a comparison, numbered from `first`,
# and returns the number of the next table.
.print_comparison <- function(x, digits, first = 1L) {
  fmt <- function(v, d = digits) .format_numeric(v, d)
  label <- function(tab) .short_names(rownames(tab), 10L)
  n_table <- first
  # This function prints one table with its note and moves to the next number.
  show <- function(out, title, note) {
    cat("\n")
    width <- .print_styled_table(out, title = sprintf("Table %d. %s", n_table, title))
    .print_note(note, width)
    n_table <<- n_table + 1L
  }

  if (!is.null(x$fit)) {
    tab <- x$fit
    out <- data.frame(Model = label(tab), lapply(as.data.frame(tab), fmt, d = 3),
                      check.names = FALSE)
    show(out, "Bayesian Fit Measures",
         c("Posterior means, ordered by BRMSEA. Chi2_ppp: posterior",
           "predictive p-value (ideal > .05)."))
  }

  if (!is.null(x$loo)) {
    tab <- x$loo$table
    best <- !is.na(tab[, "elpd_diff"]) & tab[, "elpd_diff"] == 0
    out <- data.frame(Model = label(tab), ELPD = fmt(tab[, "ELPD"]), SE = fmt(tab[, "SE"]),
                      elpd_diff = fmt(tab[, "elpd_diff"]),
                      se_diff = ifelse(best, "", fmt(tab[, "se_diff"])),
                      p_worse = ifelse(best, "", fmt(tab[, "p_worse"], 3)),
                      w_stack = fmt(tab[, "w_stacking"], 3),
                      w_pbma = fmt(tab[, "w_pseudobma"], 3))
    show(out, "LOO Comparison",
         c("Ordered by ELPD, with differences from the best model; p_worse is",
           "the probability of predicting worse than it. Weights: stacking and",
           if (x$pseudobma_bb) "pseudo-BMA+." else "pseudo-BMA."))
  }

  if (!is.null(x$ic)) {
    tab <- x$ic
    is_weight <- grepl("^w_", colnames(tab))
    cols <- lapply(seq_len(ncol(tab)), function(j) fmt(tab[, j], if (is_weight[j]) 3 else digits))
    out <- data.frame(Model = label(tab), stats::setNames(cols, colnames(tab)),
                      check.names = FALSE)
    show(out, "Information Criteria",
         "Ordered by WAIC; lower values are better. w: Akaike weights.")
  }

  if (!is.null(x$evidence$table)) {
    tab <- x$evidence$table
    out <- data.frame(Model = label(tab), log_ML = fmt(tab[, "log_ML"]),
                      MC_SE = fmt(tab[, "MC_SE"], 3), `P(M|y)` = fmt(tab[, "Estimate"], 3),
                      CI_Low = fmt(tab[, "CI_Low"], 3), CI_High = fmt(tab[, "CI_High"], 3),
                      check.names = FALSE)
    equal <- length(unique(round(x$evidence$model_prior, 12))) == 1L
    show(out, "Model Evidence",
         c("Ordered by log ML (bridge sampling). MC_SE: Micaletto and Vehtari",
           "(2025), reliable below 0.3.",
           if (equal) "Equal prior model probabilities;" else
             "Prior model probabilities as given;",
           "intervals propagate the Monte Carlo error."))
  }

  for (note in x$notes) {
    cat("\n")
    .print_note(note, 60L)
  }
  n_table
}

# This function returns the names of the models: those given, and
# Model k for the missing ones, made unique.
.model_names <- function(names, n) {
  default <- paste("Model", seq_len(n))
  if (is.null(names)) return(default)
  missing <- is.na(names) | !nzchar(names)
  names[missing] <- default[missing]
  make.unique(names, sep = "_")
}

# This function checks that the models were fitted to the same data, with
# the same data model: their log-likelihoods are comparable only then.
.check_same_data <- function(models) {
  ref <- models[[1L]]
  for (m in models[-1L]) {
    if (m$model_type != ref$model_type) {
      .befa_stop("The models must share the data model (cor, cov or raw); ",
                 "their log-likelihoods describe different data otherwise.")
    }
    same_values <- if (ref$data$summary) {
      isTRUE(all.equal(m$data$S, ref$data$S)) &&
        isTRUE(all.equal(m$data$y_bar, ref$data$y_bar))
    } else {
      isTRUE(all.equal(m$data$Y, ref$data$Y))
    }
    same <- identical(m$data$summary, ref$data$summary) &&
      identical(m$data$N, ref$data$N) && identical(m$data$J, ref$data$J) &&
      same_values
    if (!same) {
      .befa_stop("The models must be fitted to the same data, with the same ",
                 "rows and items.")
    }
  }
  invisible(TRUE)
}

# This function checks the prior model probabilities and normalizes them;
# NULL gives equal ones.
.model_prior <- function(prior, K) {
  if (is.null(prior)) return(rep(1 / K, K))
  if (!is.numeric(prior) || length(prior) != K || any(!is.finite(prior)) ||
      any(prior <= 0)) {
    .befa_stop("`model_prior` must hold ", K, " positive numbers, one per model.")
  }
  prior / sum(prior)
}

# This function returns the fit measures of a fit (when `fit` asks for the
# fit indices), its loo object and its information criteria: those stored in
# its fit indices, or computed from its draws (without the null model when
# the fit indices are not needed).
.compare_measures <- function(object, fit = FALSE) {
  fm <- object$fit_indices
  if (is.null(fm) && fit) fm <- befa_fit_measures(object)
  if (!is.null(fm)) {
    return(list(fit = fm, loo = fm$loo_object, information = fm$information))
  }
  metrics <- .posterior_metrics(object)
  loo_config <- object$options$loo_args
  if (is.null(loo_config)) loo_config <- validate_loo_args(NULL)
  loo_res <- .loo_fit(metrics$log_lik, object$method, loo_config)
  list(loo = loo_res,
       information = .information_criteria(object, metrics$log_lik, loo_res))
}

# This function returns the posterior means of the fit indices of every
# model, with their posterior predictive p-value, ordered by BRMSEA.
.compare_fit <- function(measures) {
  vars <- c("BRMSEA", "BGamma", "Adj_BGamma", "BMc", "SRMR", "BCFI", "BTLI")
  tab <- t(vapply(measures, function(m) {
    c(Chi2_ppp = m$fit$details$chi2_ppp,
      colMeans(.stacked_draws(m$fit$fit_draws, vars)))
  }, numeric(length(vars) + 1L)))
  colnames(tab) <- c("Chi2_ppp", vars)
  tab[order(tab[, "BRMSEA"]), , drop = FALSE]
}

# This function returns the LOO comparison of the models, ordered by ELPD:
# each ELPD with its SE, its difference from the best model with its SE and
# p_worse, and the stacking and pseudo-BMA weights.
.compare_loo <- function(loos, pseudobma_bb) {
  # loo_compare() names the models in its row names or, in recent versions
  # of loo, in a `model` column.
  cmp <- as.data.frame(loo::loo_compare(loos))
  ids <- if ("model" %in% names(cmp)) as.character(cmp$model) else rownames(cmp)
  row <- match(names(loos), ids)
  stacking <- suppressWarnings(loo::loo_model_weights(loos, method = "stacking"))
  pbma <- suppressWarnings(loo::loo_model_weights(loos, method = "pseudobma",
                                                  BB = pseudobma_bb))
  elpd <- t(vapply(loos, function(l) l$estimates["elpd_loo", ], numeric(2)))
  diff <- cmp$elpd_diff[row]
  se <- cmp$se_diff[row]
  table <- cbind(ELPD = elpd[, 1L], SE = elpd[, 2L], elpd_diff = diff,
                 se_diff = se,
                 p_worse = ifelse(diff == 0, NA_real_, stats::pnorm(0, diff, se)),
                 w_stacking = as.numeric(stacking), w_pseudobma = as.numeric(pbma))
  rownames(table) <- names(loos)
  list(loos = loos, table = table[order(-table[, "ELPD"]), , drop = FALSE])
}

# This function returns the information criteria of every model with their
# Akaike weights, ordered by WAIC.
.compare_ic <- function(measures) {
  rows <- c("WAIC", "DIC", "DIC_i", "BPIC")
  ic <- t(vapply(measures, function(m) m$information[rows, "Estimate"],
                 numeric(length(rows))))
  colnames(ic) <- rows
  weights <- apply(ic, 2L, .akaike_weights)
  tab <- do.call(cbind, lapply(rows, function(r) {
    stats::setNames(cbind(ic[, r], weights[, r]), c(r, paste0("w_", r)))
  }))
  colnames(tab) <- as.vector(rbind(rows, paste0("w_", rows)))
  rownames(tab) <- rownames(ic)
  tab[order(tab[, "WAIC"]), , drop = FALSE]
}

# This function returns the Akaike weights of the values of a criterion
# (lower is better): exp(-delta / 2), normalized.
.akaike_weights <- function(x) {
  w <- exp(-0.5 * (x - min(x)))
  w / sum(w)
}

# This function runs bridge sampling on every model, and returns the bridge
# objects, the log evidence of each model and its Monte Carlo standard
# error (Micaletto & Vehtari, 2025); with repetitions, their medians. A
# model whose bridge sampling fails is noted and has no evidence.
.compare_evidence <- function(models, repetitions, bridge_args) {
  if (requireNamespace("bridgesampling", quietly = TRUE) &&
      utils::packageVersion("bridgesampling") < "1.2.1") {
    .befa_stop("The Monte Carlo error of the evidence needs bridgesampling ",
               "1.2.1 or later; update it with install.packages(\"bridgesampling\").")
  }
  bridges <- list()
  notes <- character(0)
  for (k in names(models)) {
    args <- c(list(models[[k]], repetitions = repetitions), bridge_args)
    bridges[[k]] <- tryCatch(do.call(befa_bridge_sampler, args), error = function(e) {
      notes <<- c(notes, paste0("No evidence for ", k, ": ", conditionMessage(e)))
      NULL
    })
  }
  log_ml <- vapply(names(models), function(k) {
    if (is.null(bridges[[k]])) NA_real_ else stats::median(bridges[[k]]$logml)
  }, 0)
  mc_se <- vapply(names(models), function(k) {
    mcse <- bridges[[k]]$mcse_logml
    if (is.null(mcse)) NA_real_ else stats::median(mcse)
  }, 0)
  unreliable <- names(models)[!is.na(mc_se) & mc_se > 0.3]
  if (length(unreliable)) {
    notes <- c(notes, paste0("The Monte Carlo error of the evidence of ",
                             paste(unreliable, collapse = ", "), " is above 0.3, ",
                             "so it is not reliable; run more draws."))
  }
  list(bridges = bridges, log_ml = log_ml, mc_se = mc_se, notes = notes)
}

# This function returns the posterior model probabilities with an interval
# that propagates the Monte Carlo error of the log evidences, simulated as
# normal around their estimates.
.posterior_model_probs <- function(log_ml, mc_se, prior, nsim = 4000L) {
  softmax <- function(x) {
    x <- x - max(x)
    exp(x) / sum(exp(x))
  }
  estimate <- softmax(log_ml + log(prior))
  if (any(!is.finite(mc_se))) {
    return(cbind(Estimate = estimate, CI_Low = NA_real_, CI_High = NA_real_))
  }
  sims <- vapply(seq_len(nsim), function(i) {
    softmax(stats::rnorm(length(log_ml), log_ml, mc_se) + log(prior))
  }, numeric(length(log_ml)))
  sims <- matrix(sims, nrow = length(log_ml))
  q <- apply(sims, 1L, stats::quantile, probs = c(0.025, 0.975), names = FALSE)
  cbind(Estimate = estimate, CI_Low = q[1L, ], CI_High = q[2L, ])
}
