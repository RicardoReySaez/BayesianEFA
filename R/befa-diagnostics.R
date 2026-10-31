# befa-diagnostics.R
# Check the convergence of a fit on the parameters it returns
# Author: Ricardo Rey-Sáez
# Last modified: 31-10-2026

# The checks use the final parameters only: the aligned loadings, the
# uniquenesses, the communalities and the scale and mean parameters. The Stan
# parameters behind them (the rows of Z, or the loadings before alignment)
# are not identified, and their R-hat says nothing about the fit.

# Internal functions

# This function collects the diagnostics of a fit. For MCMC: the divergent
# transitions, the transitions that reached the maximum tree depth, the
# largest R-hat and the smallest bulk and tail effective sample sizes. For
# ADVI and Pathfinder, whose draws come from an approximation: the Pareto k
# of the importance ratios between the target and the approximation, when the
# backend reports both densities. Each problem is named in `issues`.
.befa_diagnostics <- function(draws, fit, settings) {
  out <- list(method = fit$method, backend = fit$backend, max_rhat = NA_real_,
              min_ess_bulk = NA_real_, min_ess_tail = NA_real_,
              divergent = NA_integer_, treedepth = NA_integer_,
              pareto_k = NA_real_)
  if (fit$method == "mcmc") {
    summary <- posterior::summarise_draws(draws, "rhat", "ess_bulk", "ess_tail")
    out$max_rhat <- max(summary$rhat, na.rm = TRUE)
    out$min_ess_bulk <- min(summary$ess_bulk, na.rm = TRUE)
    out$min_ess_tail <- min(summary$ess_tail, na.rm = TRUE)
    if (!is.null(fit$sampler)) {
      divergent <- .stacked_draws(fit$sampler, "divergent__")
      depth <- .stacked_draws(fit$sampler, "treedepth__")
      out$divergent <- as.integer(sum(divergent))
      out$treedepth <- as.integer(sum(depth >= settings$max_treedepth))
    }
  } else {
    out$pareto_k <- .approximation_k(fit$draws)
  }
  out$issues <- .diagnostic_issues(out)
  out
}

# This function returns the Pareto k of the importance ratios between the
# target density (lp__) and the approximation (lp_approx__), or NA when the
# backend does not report both.
.approximation_k <- function(draws) {
  vars <- posterior::variables(draws)
  if (!all(c("lp__", "lp_approx__") %in% vars)) return(NA_real_)
  ratios <- as.vector(.stacked_draws(draws, "lp__") -
                        .stacked_draws(draws, "lp_approx__"))
  if (!all(is.finite(ratios))) return(NA_real_)
  psis <- suppressWarnings(loo::psis(ratios, r_eff = 1))
  as.numeric(psis$diagnostics$pareto_k)
}

# This function names the problems of a fit, with the thresholds of Vehtari
# et al. (2021) for R-hat (1.01) and the effective sample size (400), and of
# Vehtari et al. (2024) for the Pareto k (0.7).
.diagnostic_issues <- function(d) {
  c(
    if (isTRUE(d$divergent > 0))
      sprintf("%d divergent transition%s", d$divergent, if (d$divergent > 1) "s" else ""),
    if (isTRUE(d$treedepth > 0))
      sprintf("%d transition%s hit the maximum tree depth", d$treedepth,
              if (d$treedepth > 1) "s" else ""),
    if (isTRUE(d$max_rhat > 1.01)) sprintf("R-hat reaches %.3f", d$max_rhat),
    if (isTRUE(min(d$min_ess_bulk, d$min_ess_tail) < 400))
      sprintf("the effective sample size falls to %.0f",
              min(d$min_ess_bulk, d$min_ess_tail)),
    if (isTRUE(d$pareto_k > 0.7))
      sprintf("the Pareto k of the approximation is %.2f", d$pareto_k)
  )
}

# This function reports, after the progress bar, the checks that
# .befa_warn_issues() also turns into a warning.
.befa_report <- function(d, settings) {
  if (d$method == "mcmc") {
    n <- d$divergent
    if (isTRUE(n == 0)) {
      cli::cli_alert_success("No divergent transitions.")
    } else if (!is.na(n)) {
      cli::cli_alert_warning("{n} divergent transition{?s}.")
    }
    n <- d$treedepth
    depth <- settings$max_treedepth
    if (isTRUE(n == 0)) {
      cli::cli_alert_success("No transition reached the maximum tree depth ({depth}).")
    } else if (!is.na(n)) {
      cli::cli_alert_warning("{n} transition{?s} reached the maximum tree depth ({depth}).")
    }
    rhat <- sprintf("%.3f", d$max_rhat)
    if (d$max_rhat <= 1.01) {
      cli::cli_alert_success("Largest R-hat {rhat}, below 1.01.")
    } else {
      cli::cli_alert_warning("Largest R-hat {rhat}, above 1.01.")
    }
    smallest <- min(d$min_ess_bulk, d$min_ess_tail)
    ess <- round(smallest)
    if (smallest >= 400) {
      cli::cli_alert_success("Smallest effective sample size {ess}, above 400.")
    } else {
      cli::cli_alert_warning("Smallest effective sample size {ess}, below 400.")
    }
    return(invisible(d))
  }

  k <- d$pareto_k
  if (is.na(k)) {
    cli::cli_alert_info(paste("The {d$backend} backend does not report the",
                              "density of the approximation, so its Pareto k",
                              "is not available."))
  } else if (k <= 0.7) {
    cli::cli_alert_success("Pareto k of the approximation {sprintf('%.2f', k)}, below 0.7.")
  } else {
    cli::cli_alert_warning("Pareto k of the approximation {sprintf('%.2f', k)}, above 0.7.")
  }
  invisible(d)
}

# This function repeats the problems of a fit as a warning, so code that runs
# quietly still sees them, and gives the remedy of each one.
.befa_warn_issues <- function(d) {
  if (!length(d$issues)) return(invisible(FALSE))
  remedies <- c(
    if (isTRUE(d$divergent > 0)) "raise `adapt_delta` towards 1",
    if (isTRUE(d$treedepth > 0)) "raise `max_treedepth`",
    if (isTRUE(d$max_rhat > 1.01) ||
        isTRUE(min(d$min_ess_bulk, d$min_ess_tail) < 400))
      "increase `iter_sampling`",
    if (isTRUE(d$pareto_k > 0.7)) "use `method = \"mcmc\"`"
  )
  .befa_warn("the posterior may be unreliable: ",
             paste(d$issues, collapse = "; "), ".\n",
             "  To fix it, ", paste(remedies, collapse = ", or "), ".")
  invisible(TRUE)
}
