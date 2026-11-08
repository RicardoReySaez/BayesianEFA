# befa-list.R
# Fit, compare and report one model per number of factors
# Author: Ricardo Rey-Sáez
# Last modified: 08-11-2026

# Public functions

#' Print a List of BEFA Models
#'
#' Prints the overview of the models fitted by [befa()] with several values
#' of `n_factors`, one row per model.
#'
#' @param x A `befa_list` object.
#' @param digits Number of decimal places of the ELPD and the WAIC.
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#' @method print befa_list
#' @export
print.befa_list <- function(x, digits = 2, ...) {
  .print_list_header(.list_header(x))
  overview <- .list_overview(x, attr(x, "comparison"))
  width <- .print_overview(overview, digits, "Overview of the Models")
  .print_note(c("summary() shows the comparison and the loadings;",
                sprintf("x$%s is the fit of that model.", names(x)[1L])), width)
  .print_failed(attr(x, "failed"), width)
  invisible(x)
}

#' Summarize a List of BEFA Models
#'
#' Summarises the models fitted by [befa()] with several values of
#' `n_factors`: an overview with one row per model, the tables of
#' [befa_compare()] when the models were compared, and the loadings and the
#' explained variance of one of them.
#'
#' @param object A `befa_list` object.
#' @param show The models whose loadings are shown: `"best"` (default), the
#'   one with the highest ELPD (without an ELPD, the highest posterior model
#'   probability); `"all"`; or some of the numbers of factors fitted, such as
#'   `2` or `c(2, 3)`. Showing a model is not selecting it: the comparison
#'   tables hold the evidence for that decision.
#' @param probs,cutoff,sort,signif_stars,name_width See [summary.befa()].
#' @param ... Ignored.
#'
#' @return A `summary.befa_list` object with `$header`, `$overview` (one row
#'   per model), `$comparison`, `$shown` (the models whose loadings are
#'   shown), `$models` (their [summary.befa()]), and `$best`, `$simpler`,
#'   `$issues` and `$failed`, which [print.summary.befa_list()] reports.
#'
#' @examples
#' \dontrun{
#' fits <- befa(HS_data, n_factors = 1:3, seed = 17)
#' summary(fits)
#' summary(fits, show = 2, cutoff = 0.3)
#' }
#'
#' @method summary befa_list
#' @export
summary.befa_list <- function(object, show = "best", probs = c(0.025, 0.975),
                              cutoff = 0, sort = TRUE, signif_stars = FALSE,
                              name_width = 8L, ...) {
  comparison <- attr(object, "comparison")
  overview <- .list_overview(object, comparison)
  best <- .best_model(overview)
  shown <- .shown_models(show, overview$Model, best$model)
  models <- lapply(object[shown], summary, probs = probs, cutoff = cutoff,
                   sort = sort, signif_stars = signif_stars,
                   name_width = name_width)
  issues <- unlist(lapply(names(object), function(k) {
    found <- object[[k]]$diagnostics$issues
    if (length(found)) paste0(k, ": ", paste(found, collapse = "; "))
  }))
  structure(list(
    header = .list_header(object), overview = overview, comparison = comparison,
    shown = shown, models = models,
    best = if (identical(show, "best")) best,
    simpler = if (identical(show, "best")) .simpler_models(comparison, best$model),
    issues = issues, failed = attr(object, "failed")
  ), class = "summary.befa_list")
}

#' Print the Summary of a List of BEFA Models
#'
#' Prints the tables of [summary.befa_list()].
#'
#' @param x A `summary.befa_list` object.
#' @param digits Number of decimal places.
#' @param cutoff,sort,signif_stars,name_width Override the options given to
#'   [summary.befa_list()].
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#' @method print summary.befa_list
#' @export
print.summary.befa_list <- function(x, digits = 2, cutoff = NULL, sort = NULL,
                                    signif_stars = NULL, name_width = NULL, ...) {
  .print_list_header(x$header)
  title <- "Table 1. Overview of the Models"
  width <- .print_overview(x$overview, digits, title)
  n_table <- 2L
  if (!is.null(x$comparison)) {
    n_table <- .print_comparison(x$comparison, digits, first = n_table)
  }

  if (!is.null(x$best)) {
    cat("\n")
    .print_note(.shown_note(x$best, x$simpler), width)
  }
  for (k in x$shown) {
    options <- .loading_options(x$models[[k]], cutoff, sort, signif_stars,
                                name_width)
    n_table <- .print_loadings(x$models[[k]], digits, options, n_table,
                               label = sprintf(" (%s)", k))$next_table
  }

  if (length(x$issues)) {
    cat("\n")
    .print_note(c("the posterior may be unreliable in",
                  paste(x$issues, collapse = "; ")), width, prefix = "Warning: ")
  }
  .print_failed(x$failed, width)
  invisible(x)
}

# Internal functions

# This function warns once that a quantity cannot use summary statistics.
.warn_no_raw_data <- function(what) {
  .befa_warn(what, " need the raw `data`; they are skipped with summary ",
             "statistics.")
}

# This function fits one model per number of factors with the same settings
# and seed, sharing the saturated and null models, and compares them. A model
# that fails is reported and left out, and the others go on. With summary
# statistics, what needs the raw data is skipped with a single warning.
.befa_fit_list <- function(inputs, run, criteria) {
  call <- run$call
  if (run$prepared$summary) {
    if (run$factor_scores) .warn_no_raw_data("Factor scores")
    if (run$compute_fit_indices) .warn_no_raw_data("Fit indices")
    run$factor_scores <- FALSE
    run$compute_fit_indices <- FALSE
    criteria <- intersect(criteria, "evidence")
  }
  if (!run$compute_fit_indices) criteria <- setdiff(criteria, "fit")

  shared <- new.env(parent = emptyenv())
  fits <- list()
  failed <- character(0)
  for (k in seq_along(inputs)) {
    input <- inputs[[k]]
    M <- input$n_factors
    label <- paste0("nf", M)
    run$call$n_factors <- M
    # The title heads the progress bar, which `refresh = 0` hides.
    run$title <- sprintf("Estimating model %d of %d (%d factor%s)", k,
                         length(inputs), M, if (M == 1L) "" else "s")
    if (run$verbose && k > 1L) cli::cli_text("")
    if (run$verbose && run$settings$refresh == 0) cli::cli_alert_info(run$title)
    fit <- tryCatch(.with_model_label(.befa_fit_one(input, run, shared), label),
                    error = function(e) conditionMessage(e))
    if (is.character(fit)) failed[label] <- fit else fits[[label]] <- fit
  }
  if (!length(fits)) {
    .befa_stop("No model could be fitted.\n",
               paste0("  ", names(failed), ": ", failed, collapse = "\n"))
  }
  if (length(failed)) {
    .befa_warn("Models not fitted: ",
               paste0(names(failed), " (", failed, ")", collapse = "; "), ".")
  }

  comparison <- NULL
  if (length(fits) > 1L && length(criteria)) {
    if (run$verbose) cli::cli_alert_info("Comparing the models.")
    comparison <- tryCatch(
      befa_compare(fits, criteria = criteria, seed = run$settings$seed),
      error = function(e) {
        .befa_warn("The models could not be compared: ", conditionMessage(e))
        NULL
      })
  }
  structure(fits, class = c("befa_list", "list"),
            n_factors = vapply(fits, `[[`, 1L, "n_factors", USE.NAMES = FALSE),
            comparison = comparison, null_model = shared$null, call = call,
            failed = failed)
}

# This function evaluates the fit of a model, naming that model in the
# warnings it raises.
.with_model_label <- function(expr, label) {
  withCallingHandlers(expr, warning = function(w) {
    .befa_warn(label, ": ", conditionMessage(w))
    invokeRestart("muffleWarning")
  })
}

# This function returns what the models of a list have in common, for the
# header of its tables.
.list_header <- function(fits) {
  first <- fits[[1L]]
  list(n_factors = attr(fits, "n_factors"), model_type = first$model_type,
       lambda_prior = first$lambda_prior, N = first$data$N,
       rotation = first$rotation, method = first$method,
       backend = first$backend)
}

# This function prints the header of the tables of a list of models.
.print_list_header <- function(h) {
  cat(sprintf("\nBayesian EFA with %s factors\n",
              paste(h$n_factors, collapse = ", ")))
  .print_note(sprintf("Data model: %s; prior: %s; N = %d; rotation: %s; %s with %s.",
                      h$model_type, h$lambda_prior, h$N, h$rotation,
                      .method_name(h$method), h$backend), 80L, prefix = "")
}

# This function returns one row per model with its convergence and its main
# indicators: the posterior means of BRMSEA and BCFI, the ELPD, the WAIC and,
# when the evidence was compared, the posterior probability of the model.
# Columns with no value are left out.
.list_overview <- function(fits, comparison) {
  diagnostic <- function(name) {
    vapply(fits, function(f) as.numeric(f$diagnostics[[name]]), 0)
  }
  index <- function(name) {
    vapply(fits, function(f) {
      fm <- f$fit_indices
      if (is.null(fm)) NA_real_ else mean(fm$posterior_fit[[name]])
    }, 0)
  }
  criterion <- function(name, table) {
    if (!is.null(table)) return(unname(table[names(fits), name]))
    vapply(fits, function(f) {
      ic <- f$fit_indices$information
      if (is.null(ic)) NA_real_ else ic[name, "Estimate"]
    }, 0)
  }
  evidence <- comparison$evidence$table
  out <- data.frame(
    Model = names(fits), Rhat = diagnostic("max_rhat"),
    ESS = pmin(diagnostic("min_ess_bulk"), diagnostic("min_ess_tail")),
    Div = diagnostic("divergent"), `Pareto k` = diagnostic("pareto_k"),
    BRMSEA = index("BRMSEA"), BCFI = index("BCFI"),
    ELPD = criterion("ELPD", comparison$loo$table),
    WAIC = criterion("WAIC", comparison$ic),
    `P(M|y)` = if (is.null(evidence)) NA_real_ else unname(evidence[names(fits), "Estimate"]),
    check.names = FALSE, row.names = NULL
  )
  out[, c(TRUE, vapply(out[-1L], function(v) any(!is.na(v)), TRUE)), drop = FALSE]
}

# This function prints the overview of the models and returns its width.
.print_overview <- function(x, digits, title) {
  decimals <- c(Rhat = 3, ESS = 0, Div = 0, `Pareto k` = 2, BRMSEA = 3, BCFI = 3,
                ELPD = digits, WAIC = digits, `P(M|y)` = 3)
  out <- x
  for (column in names(x)[-1L]) {
    out[[column]] <- .format_numeric(x[[column]], decimals[[column]])
  }
  cat("\n")
  width <- .print_styled_table(out, title = title)
  has <- function(column) column %in% names(x)
  .print_note(c(
    if (has("Rhat")) "Rhat and ESS: worst values across the parameters.",
    if (has("Div")) "Div: divergent transitions.",
    if (has("Pareto k")) "Pareto k of the approximation (ideal < 0.7).",
    if (has("BRMSEA")) "BRMSEA and BCFI: posterior means."
  ), width)
  invisible(width)
}

# This function returns the model with the highest ELPD and the name of that
# criterion, or the highest posterior probability when there is no ELPD.
.best_model <- function(overview) {
  for (column in c("ELPD", "P(M|y)")) {
    if (column %in% names(overview)) {
      return(list(model = overview$Model[which.max(overview[[column]])],
                  criterion = column))
    }
  }
  list(model = character(0), criterion = NA_character_)
}

# This function returns the names of the models to show, in the order they
# were fitted, from the `show` argument of summary.befa_list().
.shown_models <- function(show, models, best) {
  if (identical(show, "best")) return(best)
  if (identical(show, "all")) return(models)
  chosen <- if (is.numeric(show)) paste0("nf", show) else show
  if (!is.character(chosen) || !length(chosen) ||
      length(setdiff(chosen, models))) {
    .befa_stop("`show` must be \"best\", \"all\" or numbers of factors among ",
               "those fitted (", paste(sub("^nf", "", models), collapse = ", "),
               ").")
  }
  intersect(models, chosen)
}

# This function returns the models with fewer factors than `model` whose ELPD
# is within two standard errors of the difference from the best one.
.simpler_models <- function(comparison, model) {
  tab <- comparison$loo$table
  if (is.null(tab) || !length(model)) return(character(0))
  factors <- function(x) as.integer(sub("^nf", "", x))
  close <- tab[, "elpd_diff"] + 2 * tab[, "se_diff"] >= 0
  simpler <- rownames(tab)[close & factors(rownames(tab)) < factors(model)]
  simpler[order(factors(simpler))]
}

# This function explains which model the summary shows and how to see others.
.shown_note <- function(best, simpler) {
  if (!length(best$model)) {
    return(paste("No model has an ELPD to choose from; `show` chooses one, as",
                 "in summary(x, show = 2)."))
  }
  label <- if (best$criterion == "ELPD") "ELPD" else "posterior probability"
  c(sprintf("Only %s is shown, the model with the highest %s.", best$model, label),
    if (length(simpler)) {
      sprintf("%s fewer factors and an ELPD within 2 SE of it.",
              paste0(paste(simpler, collapse = ", "),
                     if (length(simpler) > 1L) " have" else " has"))
    },
    "`show` chooses others: summary(x, show = \"all\").")
}

# This function prints the models that could not be fitted, if any.
.print_failed <- function(failed, width) {
  if (!length(failed)) return(invisible(NULL))
  cat("\n")
  .print_note(paste0(names(failed), " (", failed, ")", collapse = "; "), width,
              prefix = "Not fitted: ")
}
