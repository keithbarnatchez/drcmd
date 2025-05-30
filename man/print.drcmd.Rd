% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{print.drcmd}
\alias{print.drcmd}
\title{Print drcmd object}
\usage{
\method{print}{drcmd}(x, ...)
}
\arguments{
\item{x}{An object of class drcmd}
}
\description{
S3 method for printing drcmd objects. Provides concise summary
of results, and prints values of optional arguments. Use summary() function for
more detailed summary of results.
}
\examples{
\dontrun{
n <- 2500
X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
Y <-  rbinom(n,1,plogis(X-A)) # rnorm(n) + A + X + X^2 + A*X + sin(X) # note: true ATE is 1
Ystar <- Y + rnorm(n)/2 ; R <- rbinom(n,1,plogis(X)) # error-prone outcome measurements

# Make A NA if R==0
A[R==0] <- NA
covariates <- data.frame(X1=X,X2=X2)

# Obtain ATE estimates, fitting all nuisance models with ensemble of splines +
# GAMs (save for the pseudo-outcome regression, which is done with XGboost)
drcmd_res <- drcmd(Y,A,covariates,
                   default_learners= c('SL.gam','SL.glm'),
                   po_learners = 'SL.gam')
print(drcmd_res)
}
}
