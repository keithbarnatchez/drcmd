# drcmd_demo.R
# space for testing functions in the drcmd package
#-------------------------------------------------------------------------------
library(foreach)
library(doParallel)
library(SuperLearner)
#-------------------------------------------------------------------------------
# Parameters for drcmd

default_learners <- c('SL.glm', 'SL.gam') # fit nuisance functions with ensemble of GLMS + GAMS
#-------------------------------------------------------------------------------
# Simulate a missing outcome with error-prone measurements example

n <- 1000 # sample size
X <- rnorm(n) # covariate
A <- rbinom(n,1,plogis(X)) # logit(P(A=1|X)) = X

Y <- A + X + rnorm(n) # outcome (true ATE=1)
Ystar <- Y + rnorm(n)/2 # error prone outcome measurement

R <- rbinom(n,1,0.5*plogis(X)) # missingness indicator

# Make Y NA if R==0
Y[R==0] <- NA
covariates <- data.frame(X=X)
W <- data.frame(Ystar)

# Run drcmd
drcmd_res <- drcmd::drcmd(Y,A,
               covariates,
               W=W,  
               default_learners=default_learners
               )

summary(drcmd_res)
plot(drcmd_res)
#-------------------------------------------------------------------------------
# Additional params

# Flexible choices for nuisance estimation
# E.g. use GLMs for all regressions besides pot. outcome regression, where we use GAMs
drcmd_res <- drcmd::drcmd(Y,A,
                          covariates,
                          W=W,  
                          default_learners=default_learners,
                          po_learners='SL.glm'
                          )
