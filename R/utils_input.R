#' Check if Model is Unidimensional Forced (Internal)
#'
#' Determines if the model configuration forces the unidimensional ("uni") Stan model.
#' This occurs when M=1 and using unit_vector prior.
#'
#' @param n_factors Integer. Number of latent factors.
#' @param lambda_prior Character. The lambda prior type.
#' @return Logical. TRUE if the uni model should be used.
#' @keywords internal
#' @noRd
is_uni_model <- function(n_factors, lambda_prior) {
  n_factors == 1 && lambda_prior == "unit_vector"
}

#' Validate RSP Arguments (Internal)
#'
#' Validates and fills defaults for the Rotation-Sign-Permutation algorithm.
#'
#' @param rsp_args List or NULL. User provided arguments.
#' @return A named list with 'max_iter' and 'threshold'.
#' @keywords internal
#' @noRd
validate_rsp_args <- function(rsp_args) {
  # Default settings
  defaults <- list(max_iter = 1000, threshold = 1e-6)

  # If NULL, use defaults
  if (is.null(rsp_args)) {
    return(defaults)
  }

  # If user provides it, check
  if (!is.list(rsp_args)) {
    stop("Argument 'rsp_args' must be a named list (e.g., list(max_iter=1000)).", call. = FALSE)
  }

  # If names are wrong, the're just ignored
  valid_names <- names(defaults)
  unknowns <- setdiff(names(rsp_args), valid_names)
  if (length(unknowns) > 0) {
    warning(sprintf("Unknown arguments in 'rsp_args' ignored: %s", paste(unknowns, collapse = ", ")), call. = FALSE)
  }

  # Merge defaults with user input
  out <- defaults
  if (!is.null(rsp_args$max_iter)) out$max_iter <- rsp_args$max_iter
  if (!is.null(rsp_args$threshold)) out$threshold <- rsp_args$threshold

  # Checking plausible values
  if (!is.numeric(out$max_iter) || length(out$max_iter) != 1 || out$max_iter < 1) {
    stop("'rsp_args$max_iter' must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(out$threshold) || length(out$threshold) != 1 || out$threshold <= 0) {
    stop("'rsp_args$threshold' must be a positive numeric value.", call. = FALSE)
  }

  return(out)
}

#' Validate loo_args (Internal)
#'
#' Validates and sets defaults for loo configuration.
#'
#' @param loo_args List of arguments.
#' @return List with validated arguments.
#' @keywords internal
#' @noRd
validate_loo_args <- function(loo_args) {
  # Defaults
  config <- list(
    r_eff = FALSE,
    cores = 1
  )

  if (is.null(loo_args)) {
    return(config)
  }
  if (!is.list(loo_args)) stop("Argument 'loo_args' must be a named list.", call. = FALSE)

  # Validate and update
  input_names <- names(loo_args)
  valid_names <- names(config)
  unknown_args <- setdiff(input_names, valid_names)
  if (length(unknown_args) > 0) {
    warning("Unknown arguments in 'loo_args' ignored: ", paste(unknown_args, collapse = ", "), call. = FALSE)
  }

  if (!is.null(loo_args$r_eff)) {
    if (!is.logical(loo_args$r_eff) || length(loo_args$r_eff) != 1) {
      stop("loo_args$r_eff must be a single logical value.", call. = FALSE)
    }
    config$r_eff <- loo_args$r_eff
  }

  if (!is.null(loo_args$cores)) {
    if (!is.numeric(loo_args$cores) || length(loo_args$cores) != 1 || loo_args$cores < 1) {
      stop("loo_args$cores must be a positive integer.", call. = FALSE)
    }
    config$cores <- as.integer(loo_args$cores)
  }

  return(config)
}
