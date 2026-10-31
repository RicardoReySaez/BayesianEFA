# befa-factor-scores.R
# Posterior draws of the factor scores
# Author: Ricardo Rey-Sáez
# Last modified: 31-10-2026

# Public functions

#' Bayesian Factor Scores
#'
#' Draws the factor scores of every person from their conditional posterior
#' distribution, once per posterior draw of the (aligned) model parameters.
#'
#' @param object A `befa` object returned by [befa()], fitted to raw `data`.
#' @param post_summaries Logical. If `TRUE` (default), also returns posterior
#'   summaries computed with [posterior::summarise_draws()].
#'
#' @details
#' For person \eqn{i} and each posterior draw of \eqn{\mathbf{\Lambda}},
#' \eqn{\mathbf{\Psi}} and \eqn{\boldsymbol{\nu}}, the factor scores follow
#' \deqn{\boldsymbol{\eta}_i \mid \mathbf{y}_i \sim \mathcal{N}\left(
#' \mathbf{\Sigma}_{\eta} \mathbf{\Lambda}^\top \mathbf{\Psi}^{-1}
#' (\mathbf{y}_i - \boldsymbol{\nu}),\ \mathbf{\Sigma}_{\eta}\right),
#' \qquad \mathbf{\Sigma}_{\eta} = (\mathbf{I}_M + \mathbf{\Lambda}^\top
#' \mathbf{\Psi}^{-1} \mathbf{\Lambda})^{-1}}
#' (Levy & Mislevy, 2016, pp. 199--200), computed in the metric the model was
#' fitted to. With missing values, the scores condition on the observed items
#' of each person only. The rows are those the model used: the complete rows
#' under listwise deletion, and the rows with at least one observed value
#' under FIML. Set the seed with [set.seed()] for reproducible draws.
#'
#' @return A list of class `befa_factor_scores` with:
#' * `$eta_draws`: A `draws_array` with the scores `eta[i,m]` of person `i` on
#'   factor `m`.
#' * `$summary`: Their posterior summaries (if `post_summaries = TRUE`).
#'
#' @references
#' Levy, R., & Mislevy, R. J. (2016). *Bayesian psychometric modeling*.
#' Chapman & Hall/CRC.
#'
#' @examples
#' \dontrun{
#' fit <- befa(HS_data, n_factors = 3, seed = 17)
#' scores <- befa_factor_scores(fit)
#' matrix(scores$summary$mean, ncol = 3)
#' }
#'
#' @export
befa_factor_scores <- function(object, post_summaries = TRUE) {
  .check_befa_object(object)
  draws <- object$eta
  if (is.null(draws)) {
    if (object$data$summary) {
      .befa_stop("Factor scores need the raw `data`; this model was fitted to ",
                 "summary statistics.")
    }
    draws <- .factor_score_draws(object)
  }
  out <- list(eta_draws = draws)
  if (post_summaries) out$summary <- posterior::summarise_draws(draws)
  structure(out, class = "befa_factor_scores")
}

# Internal functions

# This function draws the factor scores of a fit with the C++ sampler, from
# the rows the model was fitted to.
.factor_score_draws <- function(object) {
  p <- .fitted_parameters(object)
  Y <- object$data$Y
  observed <- lapply(seq_len(nrow(Y)), function(i) {
    as.integer(which(!is.na(Y[i, ])) - 1L)
  })
  eta <- compute_factor_scores_cpp(p$lambda, p$psi, p$nu, Y, p$J, p$M,
                                   anyNA(Y), observed)
  N <- nrow(Y)
  colnames(eta) <- paste0("eta[", rep(seq_len(N), p$M), ",",
                          rep(seq_len(p$M), each = N), "]")
  .as_chains(eta, p$n_iter, p$n_chains)
}
