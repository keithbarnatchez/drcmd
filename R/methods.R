#' @title Print drcmd object
#' @description S3 method for printing drcmd objects. Provides concise summary
#' of results, and prints values of optional arguments
#' @param x An object of class drcmd
#' @export
print.drcmd <- function(x, ...) {
  cat("drcmd results\n")
  cat("-------------------------------------------------\n")
  cat("ATE estimate: ",x$results$estimates$psi_hat_ate,'\n')
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
#' @description S3 method for summarizng drcmd results. Provides detailed summary
#' of estimation output, while also providing information on missingness
#' mechanism
#'
#' @param results An object of class drcmd
#'
#' @return A summary of the results
#'  \describe{
#'    \item{Estimates}{A data frame containing the estimates}
#'  }
#'
#'
#' @export
summary.drcmd <- function(x, ...) {
  cat("=================================================\n")
  cat("            Summary of drcmd results             \n")
  cat("=================================================\n")

  # Define the header for the table
  cat(sprintf("%-25s %10s %10s\n", "Metric", "Estimate", "SE"))
  cat("-------------------------------------------------\n")

  # Print estimates and SEs with proper alignment
  cat(sprintf("%-25s %10.3f %10.3f\n",
              "ATE:", x$results$estimates$psi_hat_ate, x$results$ses$psi_hat_ate))
  cat(sprintf("%-25s %10.3f %10.3f\n",
              "E[Y(1)]:", x$results$estimates$psi_1_hat, x$results$ses$psi_1_hat))
  cat(sprintf("%-25s %10.3f %10.3f\n",
              "E[Y(0)]:", x$results$estimates$psi_0_hat, x$results$ses$psi_0_hat))

  cat("-------------------------------------------------\n")
  cat(sprintf("%-25s %s\n", "Variables with missingness (U):", paste(x$U, collapse = ", ")))
  cat(sprintf("%-25s %s\n", "Variables without missingness (Z):", paste(x$Z, collapse = ", ")))
  cat("-------------------------------------------------\n")
  cat("Validity of results requires causal assumptions to hold,\n")
  cat("as well as the assumption that U is independent of\n")
  cat("R given Z\n")
}


# summary.drcmd <- function(x, ...) {
#   cat("=================================================\n")
#   cat("            Summary of drcmd results             \n")
#   cat("=================================================\n")
#   cat("ATE estimate (SE):     ",x$results$estimates$psi_hat_ate,' (',
#       x$results$ses$psi_hat_ate,')\n',sep="")
#   cat('E[Y(1)] estimate (SE): ', x$results$estimates$psi_1_hat, ' (',
#       x$results$ses$psi_1_hat,')\n',sep="")
#   cat('E[Y(0)] estimate (SE): ', x$results$estimates$psi_0_hat, ' (',
#       x$results$ses$psi_0_hat,')\n',sep="")
#   cat("-------------------------------------------------\n")
#   cat("Variables with missingness (U): ", x$U ,"\n")
#   cat("-------------------------------------------------\n")
#   cat("Variables without missingness (Z): ", x$Z ,"\n")
#   cat("-------------------------------------------------\n")
#   cat('Validity of results requires causal assumptions to hold,\n')
#   cat('as well as the assumption that U is independent of\nR given Z')
# }


plot.drcmd <- function(x, type = "po") {

}
