#' The 'BayesianEFA' package.
#'
#' @description Bayesian Exploratory Factor Analysis with Stan. The Stan
#' program is generated for each model structure and run with stanli,
#' cmdstanr, stanr or rstan, by MCMC (NUTS), ADVI or Pathfinder. Posterior
#' draws are post-processed with Varimax rotation and the RSP alignment
#' algorithm to resolve rotational indeterminacy. Provides Bayesian fit
#' indices, reliability estimates (Omega), and factor scores with full
#' uncertainty quantification.
#'
#' @docType package
#' @name BayesianEFA-package
#' @aliases BayesianEFA
#' @useDynLib BayesianEFA, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom loo loo relative_eff
#' @importFrom posterior as_draws_array
#' @importFrom mvnfast dmvn rmvn
#' @importFrom stats complete.cases dnorm na.omit qnorm quantile rnorm sd
#' @importFrom withr local_seed
#'
#' @references
#' Stan Development Team (2025). Stan Modeling Language Users Guide and
#' Reference Manual. <https://mc-stan.org>
#'
"_PACKAGE"
