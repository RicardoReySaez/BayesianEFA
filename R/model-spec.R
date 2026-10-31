# model-spec.R
# Describe the structure of a BEFA model, which fixes its Stan code
# Author: Ricardo Rey-Sáez
# Last modified: 31-10-2026

# Internal functions

# This function returns the structure of a model: everything that changes its
# Stan code and nothing that only changes its data, so two fits with the same
# structure share one compiled program. `loadings` is "uv" for the unit-vector
# prior with two or more factors, "uni" for that prior with one factor (a
# scaled Beta on each loading), and "normal" for independent normal loadings
# in the metric of the items. `patterns` is "single" when every row has the
# same observed items and "multiple" otherwise. `engine` names how the
# likelihood is computed.
.befa_spec <- function(model_type = c("cor", "cov", "raw"),
                       lambda_prior = c("unit_vector", "normal"),
                       n_factors,
                       patterns = c("single", "multiple"),
                       engine = "cholesky") {
  model_type <- match.arg(model_type)
  lambda_prior <- match.arg(lambda_prior)
  patterns <- match.arg(patterns)
  engine <- match.arg(engine, "cholesky")
  .check_count(n_factors, "n_factors")

  if (lambda_prior == "normal" && model_type == "cor") {
    .befa_stop("The normal prior on the loadings is not available for the ",
               "correlation model; use `lambda_prior = \"unit_vector\"`.")
  }
  loadings <- if (lambda_prior == "normal") {
    "normal"
  } else if (n_factors == 1L) {
    "uni"
  } else {
    "uv"
  }
  structure(list(model_type = model_type, loadings = loadings,
                 patterns = patterns, engine = engine),
            class = "befa_spec")
}

# This function returns the structure of the null (independence) model that
# the incremental fit indices compare against: a diagonal covariance matrix,
# with free means in the raw model. Only the covariance and raw models need it
# in Stan, since the null correlation model has no free parameters.
.befa_null_spec <- function(model_type = c("cov", "raw"),
                            patterns = c("single", "multiple")) {
  structure(list(model_type = match.arg(model_type), loadings = "none",
                 patterns = match.arg(patterns), engine = "diagonal"),
            class = "befa_spec")
}

# This function names a structure, such as "cor-uv-single-cholesky". The name
# labels the generated code.
.spec_label <- function(spec) {
  paste(spec$model_type, spec$loadings, spec$patterns, spec$engine, sep = "-")
}

# This function tells whether a structure scales a standardized solution by
# the item standard deviations `sigma`: the unit-vector loadings outside the
# correlation model, and the null model.
.spec_has_scale <- function(spec) {
  spec$model_type != "cor" && spec$loadings != "normal"
}
