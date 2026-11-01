# posterior-draws.R
# Extract and summarise the posterior draws of a fitted model
# Author: Ricardo Rey-Sáez
# Last modified: 31-10-2026

# Public functions

#' Extract Posterior Draws
#'
#' Returns the posterior draws of a `befa` fit in one of the formats of the
#' posterior package, in the metric of the model or standardized.
#'
#' @param befa_fit A `befa` object returned by [befa()].
#' @param pars Character vector of parameter names, such as `"Lambda_std"`
#'   (all standardized loadings) or `"h2[3]"`. `NULL` (default) returns every
#'   parameter.
#' @param format Output format: `"matrix"` (default), `"array"`, `"df"`,
#'   `"list"` or `"rvars"`.
#' @param metric `"model"` (default) returns the parameters as the model
#'   estimates them; `"std"` returns the standardized loadings `Lambda_std`,
#'   the standardized uniquenesses `Psi_std` and the communalities `h2`; and
#'   `"raw"` the parameters in the units of the items, which the correlation
#'   model does not have.
#'
#' @return The posterior draws in the requested format.
#'
#' @examples
#' \dontrun{
#' fit <- befa(HS_data, n_factors = 3, model = "raw", seed = 17)
#' extract_posterior_draws(fit, pars = "Lambda", format = "rvars")
#' extract_posterior_draws(fit, pars = "Lambda_std", metric = "std")
#' }
#'
#' @export
extract_posterior_draws <- function(befa_fit, pars = NULL,
                                    format = c("matrix", "array", "df", "list",
                                               "rvars"),
                                    metric = c("model", "std", "raw")) {
  .check_befa_object(befa_fit)
  format <- match.arg(format)
  draws <- .metric_draws(befa_fit, match.arg(metric))
  if (!is.null(pars)) draws <- posterior::subset_draws(draws, variable = pars)
  switch(format,
    matrix = posterior::as_draws_matrix(draws),
    array = draws,
    df = posterior::as_draws_df(draws),
    list = posterior::as_draws_list(draws),
    rvars = posterior::as_draws_rvars(draws)
  )
}

#' Compute Posterior Summaries
#'
#' Summarises the posterior draws of a `befa` fit with
#' [posterior::summarise_draws()].
#'
#' @inheritParams extract_posterior_draws
#' @param ... Summary functions passed to [posterior::summarise_draws()]; by
#'   default its mean, median, SD, MAD, 5% and 95% quantiles, R-hat and
#'   effective sample sizes.
#'
#' @return A `draws_summary` data frame.
#'
#' @examples
#' \dontrun{
#' fit <- befa(HS_data, n_factors = 3, seed = 17)
#' posterior_summaries(fit, pars = "h2")
#' }
#'
#' @export
posterior_summaries <- function(befa_fit, pars = NULL,
                                metric = c("model", "std", "raw"), ...) {
  draws <- extract_posterior_draws(befa_fit, pars, format = "array",
                                   metric = match.arg(metric))
  posterior::summarise_draws(draws, ...)
}

# Internal functions

# This function checks that an object is a fit returned by befa().
.check_befa_object <- function(x) {
  if (!inherits(x, "befa")) {
    .befa_stop("`befa_fit` must be a fit returned by befa().")
  }
  invisible(TRUE)
}
