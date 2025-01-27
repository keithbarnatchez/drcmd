#'
#'
#'
#'
#'
#'
#'
#'
#'
print.drcmd <- function(x, ...) {
  cat("drcmd results\n")
  cat("-------------------------------------------------\n")
  cat("ATE estimate: ",x$estimates$psi_hat_ate,'\n')
  cat("-------------------------------------------------\n")
  cat("Variables with missingness (U): ", x$U ,"\n")
  cat("-------------------------------------------------\n")
  cat("Variables without missingness (Z): ", x$Z ,"\n")
  cat("-------------------------------------------------\n")
  cat('Validity of results requires causal assumptions to hold\n')
  cat('As well as the assumption that U is independent of R given Z\n')
  cat("Number of cross-fitting folds (k):", x$k, "\n")
  cat("Bootstrap samples (nboot):", x$nboot, "\n")
  cat("\nEstimates:\n")
}

#' @title Summarize results from drcmd
#' @description S3 method for summarizng drcmd results
#'
#' @param results An object of class drcmd
#'
#' @return A summary of the results
#'  \describe{
#'    \item{Estimates}{A data frame containing the estimates}
#'  }
#'
#'
#'
#'
#'
#'
summary.drcmd <- function(results, ...) {
  cat("=================================================\n")
  cat("            Summary of drcmd results             \n")
  cat("=================================================\n")
  cat("ATE estimate: ",results$estimates$psi_hat_ate,'\n')
  cat("-------------------------------------------------\n")
  cat("Variables with missingness (U): ", results$U ,"\n")
  cat("-------------------------------------------------\n")
  cat("Variables without missingness (Z): ", results$Z ,"\n")
  cat("-------------------------------------------------\n")
  cat('Validity of results requires causal assumptions to hold,\n')
  cat('as well as the assumption that U is independent of\nR given Z')
}
