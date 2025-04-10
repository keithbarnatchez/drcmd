# test_funcs.R
# space for testing functions in the drcmd package
#-------------------------------------------------------------------------------
rm(list=ls())
source('../R/utils.R')
source('../R/drcmd.R')
source('../R/nuis.R')
source('../R/methods.R')
# devtools::install_github('keithbarnatchez/drcmd')
#-------------------------------------------------------------------------------
library(foreach)
library(doParallel)
library(SuperLearner)
#-------------------------------------------------------------------------------
# Optional params for drcmd

eem_ind <- FALSE
default_learners <- 'SL.glm'
options(error=traceback)
#-------------------------------------------------------------------------------

n <- 3000
X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
Y <- A + X + rnorm(n)
Ystar <- Y + rnorm(n)/2 ; R <- rbinom(n,1,0.5*plogis(X)) # error-prone outcome measurements

# Make A NA if R==0
Y[R==0] <- NA
covariates <- data.frame(X1=X)

drcmd_tml <- drcmd(Y,A,covariates,
                   default_learners= c('SL.glm','SL.glm.interaction','SL.earth'),
                   eem_ind=F,tml=T,cutoff=0)
drcmd_tml$results$estimates

registerDoParallel(cores=8)
reslist <- foreach(1:100, .combine='rbind') %dopar% {

  n <- 5000
  X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
  Y <- rbinom(n,1,prob = 0.2 + 0.3*A) # A + X + rnorm(n)
  Ystar <- Y + rnorm(n)/2 ; R <- rbinom(n,1,0.5*plogis(X)) # error-prone outcome measurements

  # Make A NA if R==0
  Y[R==0] <- NA
  covariates <- data.frame(X1=X)

drcmd_tml <- drcmd(Y,A,covariates,
                   default_learners= c('SL.glm','SL.glm.interaction','SL.gam'),
                   r_learners='SL.glm',
                   po_learners='SL.gam',
                   eem_ind=F,cutoff=0)
res <- data.frame(res=drcmd_tml$results$estimates$psi_hat_ate)
}
mean(reslist$res)
sd(reslist$res)
#-------------------------------------------------------------------------------
# Make a couple functions for simulating simple missing outcome data structure

n <- 5000
X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
Y <-  rnorm(n) + A + X  # note: true ATE is 1
# Y <- rbinom(n,1,0.5)
Ystar <- Y + rnorm(n)/2 ; R <- rbinom(n,1,plogis(X)) # error-prone outcome measurements
X2=X+rnorm(n)

# Make A NA if R==0
A[R==0] <- NA
covariates <- data.frame(X1=X,X2=X2)

# Obtain ATE estimates, fitting all nuisance models with ensemble of splines +
# GAMs (save for the pseudo-outcome regression, which is done with XGboost)
drcmd_res <- drcmd(Y,A,covariates,
                   default_learners= c('SL.glm','SL.glm.interaction','SL.earth'),
                   eem_ind=F)

drcmd_tml <- drcmd(Y,A,covariates,
                   default_learners= c('SL.glm','SL.glm.interaction','SL.earth'),
                   eem_ind=F,tml=T)

drcmd_eem <- drcmd(Y,A,covariates,
                   default_learners= c('SL.glm','SL.glm.interaction','SL.earth'),
                   eem_ind=T)


summary(drcmd_res, detail=T)
summary(drcmd_tml, detail=T)
summary(drcmd_eem, detail=T)
plot(drcmd_res)

#-------------------------------------------------------------------------------
# Try with missing outcome and treatment

n <- 1e3
X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
Y <- rnorm(n) + A + X + X^2 + A*X + sin(X)
Ystar <- Y + rnorm(n)/2 ; R <- rbinom(n,1,plogis(X)) ; X <- as.data.frame(X)

# Make Y NA if R==0
Y[R==0] <- NA
A[R==0] <- NA

df <- data.frame(Y=Y,A=A,X=X,Ystar=Ystar,R=R)

# X = data.frame(cbind(X,X))
W = X[,0]

drcmd_res <- drcmd(Y,A,X, default_learners= c('SL.glm','SL.glm.interaction','SL.earth','SL.gam'),
                   eem_ind=FALSE,k=1,tml=TRUE)
#-------------------------------------------------------------------------------
# Missng outcome and proxy used

n <- 1e3
X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
Y <- rnorm(n) + A + X + X^2 + A*X + sin(X)
Ystar <- Y + rnorm(n)/2 ; R <- rbinom(n,1,plogis(X)) ; X <- as.data.frame(X)

# Make Y NA if R==0
Y[R==0] <- NA
A[1:50] <- NA
X[1:70,] <- NA

df <- data.frame(Y=Y,A=A,X=X,Ystar=Ystar,R=R)

W <- data.frame(Ystar)

drcmd_res <- drcmd(Y,A,X,
                   W=data.frame(Ystar),
                   default_learners= c('SL.glm','SL.glm.interaction',
                                       'SL.earth','SL.gam'),
                   eem_ind=TRUE,k=2)
summary(drcmd_res)
#-------------------------------------------------------------------------------
# Test sub functions

# Missing data pattern finder
V <- find_missing_pattern(Y,A,X,W)
Z <- V$Z ; R <- V$R ; X <- V$X ; Y <- V$Y ; A <- V$A
Rprobs <- NA

# Nuisance functions
idx <- 1:n
kappa_hat <- est_kappa(idx,Z,R,sl.lib)
m_hat <- est_m_a(idx,Y,A,X,R,kappa_hat,hal_ind,sl.lib)
g_hat <- est_g(idx,A,X,R,kappa_hat,hal_ind,sl.lib)
check <- get_nuisance_ests(idx,Y,A,X,Z,R,hal_ind,sl.lib,Rprobs)

# Check the phi maker
phi_hat  <- get_phi_hat(Y,A,X,R,g_hat,m_hat,kappa_hat,hal_ind,sl.lib)
phi_1_hat <- phi_hat$phi_1_hat ; phi_0_hat <- phi_hat$phi_0_hat

varphi_hat <- est_varphi(idx,R,Z,phi_1_hat,phi_0_hat,hal_ind,sl.lib)

est_psi(idx, R, Z, kappa_hat, phi_hat,varphi_hat)

drcmd_res <- drcmd(Y,A,X)
#-------------------------------------------------------------------------------
# Test coverage of drcmd


compute_CI <- function(est, se) {
  lower <- est - 1.96 * se
  upper <- est + 1.96 * se
  return(c(lower,upper))
}

trim <- function(p,cutoff=0.05) {

  # trim lower piece
  p <- ifelse(p<cutoff,cutoff,p)

  # trim upper piece
  p <- ifelse(p>(1-cutoff),1-cutoff,p)

  return(p)
}

expit <- function(o) {
  return(exp(o)/(1+exp(o)))
}
#-------------------------------------------------------------------------------
n <- 1000
p <- 3 # number of covariates
delta <- c(-0.1,-0.6,-0.9)
nu <- c(0.1,-0.1,0.1)
beta <- c(1,2,-2) # outcome model main effects
tau <- 1 # constant treatment effect

gamma <- rep(1,p) # interaction effects with tmt c(0.5,2.1,-1.2)
eta <- c(0.6,-0.2,0.8,0.1,-0.3)
rho=0.1
sl_lib <- c('SL.glm','SL.glm.interaction','SL.gam')
SL.rf = function(...) {
  SL.ranger(..., num.trees = 50)
}
#-------------------------------------------------------------------------------
# Brief sim emulating paper 2 params

nsim <- 250

registerDoParallel(cores=8)
results <- foreach(1:nsim,.combine=rbind) %dopar% {
  # covariates
  X <- matrix(runif(n*p),nrow=n,ncol=p)

  # treatment
  probs <- trim(expit(as.matrix(X)%*%delta))
  A <- rbinom(nrow(X),size=1,prob = probs)

  # outcome
  Y <- X%*%beta + tau*A + A*(X%*%gamma) + rnorm(nrow(X),mean=0,sd= rowSums(X)) # exp(rowSums(X)))

  # Y <- rbinom(n,1,0.5)
  # measurements
  Ystar <- Y + X%*%nu + rnorm(nrow(X),mean=0,sd=0.1)
  Astar <- ifelse(runif(length(A)) < 0.85, A, 1 - A)

  Rprobs <- rho*trim(expit(as.matrix(cbind(X,Astar,Ystar))%*%eta))/mean(trim(expit(as.matrix(cbind(X,Astar,Ystar))%*%eta)))
  R <- rbinom(nrow(X),size=1,prob = Rprobs)
  X <- as.data.frame(X)

  oracle_aipw <- AIPW::AIPW$new(Y = Y, A=A, W = X,Q.SL.library = 'SL.glm.interaction',
                                g.SL.library = 'SL.glm.interaction',
                                verbose=FALSE)$fit()$summary()
  oracle_ests <- oracle_aipw$estimates$RD['Estimate']

  tmle_res <- tmle::tmle(Y=Y,A=A,W=X,
  g.SL.library=sl_lib,
  Q.SL.library=sl_lib,
  obsWeights=as.double(R/Rprobs))
  tmle_ests <- tmle_res$estimates$ATE$psi

  # Make Y NA if R==0
  Y[R==0] <- NA
  A[R==0] <- NA

  ### *************
  ### *************
  drcmd_res_tml <- drcmd::drcmd(Y,A,X,
                         W=data.frame(Ystar,Astar),
                         default_learners= sl_lib,
                         Rprobs=Rprobs,
                         k=1,tml=T) ; drcmd_res_tml
  ### *************
  ### *************

  drcmd_res_tml$results$estimates

  drcmd_tml <- drcmd_res_tml$results$estimates$psi_hat_ate
  drcmd_tml_int <- compute_CI(drcmd_tml,
                              drcmd_res_tml$results$ses$psi_hat_ate)
  drcmd_tml_cov <- (drcmd_tml_int[1] <= 2.5) & (drcmd_tml_int[2] >= 2.5)

  # Estimate with eem
  drcmd_res_evm <-  drcmd::drcmd(Y,A,X,
                         W=data.frame(Ystar,Astar),
                         default_learners= sl_lib,
                         # r_learners='SL.glm',
                         Rprobs=Rprobs,
                         eem_ind=TRUE,k=1) ; drcmd_res_evm

  IF1 <- mean(drcmd_res_evm$results$nuis$m_1_hat) + mean(R/Rprobs*drcmd_res_evm$results$nuis$phi_1_hat)
  IF0 <- mean(drcmd_res_evm$results$nuis$m_0_hat) + mean(R/Rprobs*drcmd_res_evm$results$nuis$phi_0_hat)
  rw_aipws <- IF1 - IF0

  # Estimate without eem
  drcmd_res_mid <-  drcmd(Y,A,X,
                         W=data.frame(Ystar,Astar),
                         default_learners= sl_lib,
                         # r_learners='SL.glm',
                         Rprobs=Rprobs,
                         eem_ind=FALSE,k=1)

  # evm results
  drcmd_evm_ests <- drcmd_res_evm$results$estimates$psi_hat_ate
  drcmd_direct_ests <- drcmd_res_evm$results$estimates$psi_hat_ate_direct
  temp_int <- compute_CI(drcmd_res_evm$results$estimates$psi_hat_ate,
                         drcmd_res_evm$results$ses$psi_hat_ate)
  drcmd_evm_cov <- (temp_int[1] < 2.5) & (temp_int[2] > 2.5)

  # mid results
  drcmd_mid_ests <- drcmd_res_mid$results$estimates$psi_hat_ate
  temp_int <- compute_CI(drcmd_res_mid$results$estimates$psi_hat_ate,
                         drcmd_res_mid$results$ses$psi_hat_ate)
  drcmd_mid_cov <- (temp_int[1] < 2.5) & (temp_int[2] > 2.5)

  return(data.frame(oracle_ests=oracle_ests,drcmd_evm_ests=drcmd_evm_ests,
                    drcmd_mid_ests=drcmd_mid_ests,
                    drcmd_direct_ests=drcmd_direct_ests,
                    tmle_ests=tmle_ests, rw_aipws=rw_aipws,
                    drcmd_tml=drcmd_res_tml$results$estimates$psi_hat_ate,
                    drcmd_evm_cov=drcmd_evm_cov,drcmd_mid_cov=drcmd_mid_cov))
}


# Make a table with mean and SD of the above things
tbl <- data.frame(method=c('oracle','drcmd evm', 'drcmd mid', 'drcmd direct', 'drcmd tml','rw tml', 'rw aipw'),
                   mean=c(mean(results$oracle_ests),mean(results$drcmd_evm_ests),mean(results$drcmd_mid_ests),
                         mean(results$drcmd_direct_ests),mean(results$drcmd_tml),mean(results$tmle_ests),mean(results$rw_aipws)),
                   sd=c(sd(results$oracle_ests),sd(results$drcmd_evm_ests),sd(results$drcmd_mid_ests),
                        sd(results$drcmd_direct_ests),sd(results$drcmd_tml),sd(results$tmle_ests),sd(results$rw_aipws)))
tbl

mean(results$drcmd_evm_cov) ; sd(results$drcmd_evm_cov)
mean(results$drcmd_mid_cov) ; sd(results$drcmd_mid_cov)
mean(results$tmle_ests)

