# befa-reliability.R
# Omega reliability coefficients from the posterior draws
# Author: Ricardo Rey-Sáez
# Last modified: 31-10-2026

# Public functions

#' Bayesian Reliability Estimation (McDonald's Omega)
#'
#' Computes McDonald's omega for the total score and for each factor at every
#' posterior draw, so the coefficients come with their full posterior
#' distribution.
#'
#' @param object A `befa` object returned by [befa()].
#' @param probs Numeric vector of length 2 with the quantiles of the credible
#'   intervals (default `c(0.025, 0.975)`).
#'
#' @details
#' Unlike Cronbach's \eqn{\alpha}, which assumes equal loadings, omega follows
#' the congeneric structure of the factor model. The coefficients use the
#' standardized solution (\eqn{\mathbf{\Lambda}} and \eqn{\mathbf{\Psi}} of
#' standardized items, with orthogonal factors), so they describe the sum of
#' the standardized items in every data model.
#'
#' * Omega total, the share of the variance of the total score due to all the
#'   factors:
#'   \deqn{\omega_t = \frac{\sum_m (\sum_j \lambda_{jm})^2}
#'   {\sum_m (\sum_j \lambda_{jm})^2 + \sum_j \psi_j}}
#' * Omega of factor m, treating the uniquenesses as measurement error:
#'   \deqn{\omega_m = \frac{(\sum_j \lambda_{jm})^2}
#'   {(\sum_j \lambda_{jm})^2 + \sum_j \psi_j}}
#'
#' @return A list of class `befa_reliability` with:
#' * `$omega_draws`: A `draws_array` with `omega_total` and `omega_F1`, ...,
#'   `omega_FM`.
#' * `$summary`: Their posterior summaries from [posterior::summarise_draws()].
#' * `$total`: The posterior mean, SD and credible interval of omega total.
#' * `$subscales`: The same for each factor, one row per factor.
#'
#' @references
#' McDonald, R. P. (1999). *Test Theory: A Unified Treatment*. Lawrence Erlbaum
#' Associates.
#'
#' @examples
#' \dontrun{
#' fit <- befa(HS_data, n_factors = 3, seed = 17)
#' befa_reliability(fit)$subscales
#' }
#'
#' @export
befa_reliability <- function(object, probs = c(0.025, 0.975)) {
  .check_befa_object(object)
  if (!is.null(object$reliability)) return(object$reliability)
  if (!is.numeric(probs) || length(probs) != 2L || any(probs <= 0 | probs >= 1)) {
    .befa_stop("`probs` must hold two probabilities between 0 and 1.")
  }

  std <- .metric_draws(object, "std")
  J <- object$stan_data$J
  M <- object$n_factors
  lambda <- .stacked_draws(std, "Lambda_std")
  psi <- rowSums(.stacked_draws(std, "Psi_std"))

  # Sum of the loadings of each factor, draw by draw (S x M)
  sums <- vapply(seq_len(M), function(m) {
    rowSums(lambda[, (m - 1L) * J + seq_len(J), drop = FALSE])
  }, numeric(nrow(lambda)))
  sums <- matrix(sums, ncol = M)
  common <- rowSums(sums^2)
  omega <- cbind(omega_total = common / (common + psi),
                 sums^2 / (sums^2 + psi))
  colnames(omega) <- c("omega_total", paste0("omega_F", seq_len(M)))

  dims <- dim(object$draws)
  draws <- .as_chains(omega, dims[1L], dims[2L])
  describe <- function(x) {
    q <- stats::quantile(x, probs, names = FALSE)
    data.frame(Estimate = mean(x), SD = stats::sd(x), Lower = q[1L],
               Upper = q[2L])
  }
  subscales <- do.call(rbind, lapply(seq_len(M), function(m) describe(omega[, m + 1L])))
  rownames(subscales) <- paste0("F", seq_len(M))
  structure(list(
    omega_draws = draws, summary = posterior::summarise_draws(draws),
    total = as.list(describe(omega[, 1L])), subscales = subscales,
    probs = probs
  ), class = "befa_reliability")
}
