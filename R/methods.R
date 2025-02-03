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
}

#' @title Summarize results from drcmd
#' @description S3 method for summarizng drcmd results. Provides detailed summary
#' of estimation output, while also providing information on missingness
#' mechanism. Can optionally print out values of user-supplied arguments
#'
#' @param results An object of class drcmd
#'
#' @return No return value. Called for printing a detailed results summary
#'
#' @export
summary.drcmd <- function(x, detail=FALSE, ...) {
  cat("======================================================================\n")
  cat("                        Summary of drcmd results                      \n")
  cat("======================================================================\n")

  # Define the header for the table with better spacing
  cat(sprintf("%-15s %12s %12s %26s\n", "Metric", "Estimate", "SE", "95% CI"))
  cat("----------------------------------------------------------------------\n")

  # Function to compute confidence interval
  compute_CI <- function(est, se) {
    lower <- est - 1.96 * se
    upper <- est + 1.96 * se
    sprintf("[%.3f, %.3f]", lower, upper)
  }

  # Print estimates, SEs, and CIs with improved formatting
  cat(sprintf("%-15s %12.3f %12.3f %26s\n",
              "ATE:", x$results$estimates$psi_hat_ate, x$results$ses$psi_hat_ate,
              compute_CI(x$results$estimates$psi_hat_ate, x$results$ses$psi_hat_ate)))
  cat(sprintf("%-15s %12.3f %12.3f %26s\n",
              "E[Y(1)]:", x$results$estimates$psi_1_hat, x$results$ses$psi_1_hat,
              compute_CI(x$results$estimates$psi_1_hat, x$results$ses$psi_1_hat)))
  cat(sprintf("%-15s %12.3f %12.3f %26s\n",
              "E[Y(0)]:", x$results$estimates$psi_0_hat, x$results$ses$psi_0_hat,
              compute_CI(x$results$estimates$psi_0_hat, x$results$ses$psi_0_hat)))

  cat("----------------------------------------------------------------------\n")
  cat(sprintf("%-30s %s\n", "Variables with missingness (U):", paste(x$U, collapse = ", ")))
  cat(sprintf("%-30s %s\n", "Variables without missingness (Z):", paste(x$Z, collapse = ", ")))
  cat("----------------------------------------------------------------------\n")
  cat("Validity of results requires causal assumptions to hold,\n")
  cat("as well as the assumption that U is independent of R given Z\n")

  if (detail) {
    cat("----------------------------------------------------------------------\n")
    cat("Number of cross-fitting folds (k):", x$params$k, "\n")
    cat("Bootstrap samples (nboot):", x$params$nboot, "\n")
    cat('Outcome regression nuisance learners:', x$params$m_learners, '\n')
    cat('Propensity score nuisance learners:', x$params$g_learners, '\n')
    cat('Missingness nuisance learners:', x$params$r_learners, '\n')
    cat('Pseudo-outcome nuisance learners:', x$params$po_learners, '\n')
  }
}

#' @title Plot results from drcmd object
#' @description S3 method for plotting results from drcmd object. Plots are
#' available for the following: (1) psuedo outcome regression fit, (2) influence
#' curve distribution,
#'
plot.drcmd <- function(x, type = "All") {

  # First check if type is valid
  if (!(type %in% c('PO', 'IC', 'g_hat', 'r_hat','All'))) {
    stop("Invalid type argument. Must be one of 'All', PO', 'IC', 'ghat' or 'r_hat'")
  }

  if (type == "PO" | type=='All') {

    nuis <- x$results$nuis
    R <- x$R

    phis <- rbind(data.frame(phi_hat=nuis$phi_1_hat[R==1], varphi_hat=nuis$varphi_1_hat[R==1],A=1),
                  data.frame(phi_hat=nuis$phi_0_hat[R==1], varphi_hat=nuis$varphi_0_hat[R==1],A=0))


    pp <- ggplot2::ggplot(phis,
                    ggplot2::aes(x = varphi_hat, y = phi_hat-varphi_hat,group=as.factor(A),
                                 color=as.factor(A))) +
      ggplot2::geom_point(alpha=0.3) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = "Pseudo Outcome Regression Residuals vs Fitted Values",
        subtitle = 'Restricted to complete cases',
           y = "Pseudo outcome regression residuals",
           x = "Psuedo outcome regression fited values") +
      ggplot2::labs(color = "Treatment") +
    ggplot2::theme(legend.position = "bottom")
    print(pp)
    if (type=='All') {
    line <- readline(prompt = "Enter 'q' to quit or any other key to see next plot: ")
    if (line=='q') return()
    }


  }

  if (type=='IC' | type=='All') {
    IC1 <- (R*x$results$nuis$kappa_hat)*(x$results$nuis$phi_1_hat -
            x$results$nuis$varphi_1_hat) + x$results$nuis$varphi_1_hat
    ICO <- (R*x$results$nuis$kappa_hat)*(x$results$nuis$phi_0_hat -
            x$results$nuis$varphi_0_hat) + x$results$nuis$varphi_0_hat


    ICs <- rbind(data.frame(IC=IC1, A='E[Y(1)]'),
                 data.frame(IC=ICO, A='E[Y(0)]'),
                 data.frame(IC=IC1-ICO, A='ATE'))

    pp <- ggplot2::ggplot(ICs,
                    ggplot2::aes(x = IC, group=as.factor(A), color=as.factor(A))) +
      ggplot2::geom_density(alpha=0.3) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = "Influence Curve Distribution",
        y = "Density",
        x = "Influence Curve") +
      ggplot2::labs(color = "Estimand") +
      ggplot2::theme(legend.position = "bottom")
      print(pp)
      if (type=='All') {
        line <- readline(prompt = "Enter 'q' to quit or any other key to see next plot: ")
        if (line=='q') return()
      }
  }

  if (type=='g_hat' | type=='All') {
    pp <- ggplot2::ggplot(data.frame(g_hat=x$results$nuis$g_hat, A=x$R),
                    ggplot2::aes(x = g_hat, group=as.factor(A), color=as.factor(A))) +
      ggplot2::geom_density(alpha=0.3) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = "Treatment Propensity Score Distribution",
        y = "Density",
        x = "Propensity Score") +
      ggplot2::labs(color = "Treatment") +
      ggplot2::theme(legend.position = "bottom")
    print(pp)
    if (type=='All') {
      line <- readline(prompt = "Enter 'q' to quit or any other key to see next plot: ")
      if (line=='q') return()
    }
  }

  if (type=='r_hat' |  type=='All') {
    pp <- ggplot2::ggplot(data.frame(kappa_hat=x$results$nuis$kappa_hat, A=x$R),
                    ggplot2::aes(x = kappa_hat, group=as.factor(A), color=as.factor(A))) +
      ggplot2::geom_density(alpha=0.3) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = "Complete Case Probability Distribution",
        y = "Density",
        x = "Complete Case Probability") +
      ggplot2::labs(color = "Complete Case") +
      ggplot2::theme(legend.position = "bottom")
    print(pp)
  }

}







