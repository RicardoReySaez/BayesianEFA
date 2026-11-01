# precompile.R
# Knit the vignettes and the README, whose examples run Stan
# Author: Ricardo Rey-Sáez
# Last modified: 01-11-2026

# The vignettes run Stan, which CRAN machines and R CMD check may not have,
# so their sources (*.Rmd.orig) are knitted here, and the package ships the
# knitted *.Rmd with its output and figures. The README is knitted the same
# way. The package is installed into a temporary library first, so the
# chains, which run in background R processes, use this version of it.
# Run from the package root:
#   Rscript vignettes/precompile.R

lib <- file.path(tempdir(), "precompile-lib")
dir.create(lib, showWarnings = FALSE)
status <- system2(file.path(R.home("bin"), "R"),
                  c("CMD", "INSTALL", "--no-test-load",
                    paste0("--library=", shQuote(lib)), shQuote(getwd())),
                  stdout = FALSE, stderr = FALSE)
if (status != 0) stop("The package could not be installed.")

callr::r(function(root) {
  setwd(file.path(root, "vignettes"))
  for (name in c("getting_started", "alignment")) {
    knitr::knit(paste0(name, ".Rmd.orig"), paste0(name, ".Rmd"), quiet = TRUE)
  }
  setwd(root)
  rmarkdown::render("README.Rmd", quiet = TRUE)
}, args = list(root = getwd()), libpath = c(lib, .libPaths()), show = TRUE)
