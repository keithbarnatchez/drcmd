# test_funcs.R
# space for testing functions in the drcmd package
#-------------------------------------------------------------------------------
rm(list=ls())
source('../R/utils.R')
source('../R/drcmd.R')
source('../R/nuis.R')
source('../R/methods.R')
#-------------------------------------------------------------------------------
# Optional params for drcmd

eem_ind <- FALSE
default_learners <- 'hal'
#-------------------------------------------------------------------------------
# Make a couple functions for simming simple data structure

n <- 1e3
X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
Y <- rnorm(n) + A + X + X^2 + A*X + sin(X)
Ystar <- Y + rnorm(n)/2 ; R <- rbinom(n,1,plogis(X)) ; X <- as.data.frame(X)

# Make Y NA if R==0
Y[R==0] <- NA

df <- data.frame(Y=Y,A=A,X=X,Ystar=Ystar,R=R)

drcmd_res <- drcmd(Y,A,X, default_learners= c('SL.glm'),
                   po_learners='SL.hal9001',
                   eem_ind=TRUE,k=1, Rprobs=plogis(X$X))

summary(drcmd_res)
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
                   eem_ind=FALSE,k=1)
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


#-------------------------------------------------------------------------------
n <- 1000
p <- 3 # number of covariates
delta <- c(-0.1,-0.6,-0.9)
nu <- c(0.1,-0.1,0.1)
beta <- c(1,2,-2) # outcome model main effects
tau <- 1 # constant treatment effect

gamma <- rep(1,p) # interaction effects with tmt c(0.5,2.1,-1.2)
eta <- c(0.6,-0.2,0.8,0.1,-0.3)
rho=0.2
#-------------------------------------------------------------------------------
# Brief sim emulating paper 2 params

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

nsim <- 70
drcmd_evm_ests <- rep(NA,nsim)
drcmd_mid_ests <- rep(NA,nsim)
drcmd_direct_ests <- rep(NA,nsim)
oracle_ests <- rep(NA,nsim)

drcmd_evm_cov <- rep(NA,nsim)
drcmd_mid_cov <- rep(NA,nsim)
for (ss in 1:nsim) {
  print(ss)

  # covariates
  X <- matrix(runif(n*p),nrow=n,ncol=p)

  # treatment
  probs <- trim(expit(as.matrix(X)%*%delta))
  A <- rbinom(nrow(X),size=1,prob = probs)

  # outcome
  Y <- X%*%beta + tau*A + A*(X%*%gamma) + rnorm(nrow(X),mean=0,sd=1) # exp(rowSums(X)))

  # measurements
  Ystar <- Y + X%*%nu + rnorm(nrow(X),mean=0,sd=1)
  Astar <- ifelse(runif(length(A)) < 0.8, A, 1 - A)

  Rprobs <- rho*trim(expit(as.matrix(cbind(X,Astar,Ystar))%*%eta))/mean(trim(expit(as.matrix(cbind(X,Astar,Ystar))%*%eta)))
  R <- rbinom(nrow(X),size=1,prob = Rprobs)
  X <- as.data.frame(X)


  oracle_aipw <- AIPW::AIPW$new(Y = Y, A=A, W = X,Q.SL.library = 'SL.glm.interaction',
                                g.SL.library = 'SL.glm.interaction')$fit()$summary()
  oracle_ests[ss] <- oracle_aipw$estimates$RD

  # Make Y NA if R==0
  Y[R==0] <- NA
  A[R==0] <- NA

  # Estimate with eem
  drcmd_res_evm <- drcmd(Y,A,X,
                     W=data.frame(Ystar,Astar),
                     default_learners= c('SL.glm.interaction'),
                     po_learners = 'SL.gam',
                     eem_ind=TRUE,k=1,Rprobs=Rprobs) ; drcmd_res_evm

  # Estimate without eem
  drcmd_res_mid <- drcmd(Y,A,X,
                      W=data.frame(Ystar,Astar),
                      default_learners= c('SL.glm.interaction'),
                      po_learners = 'SL.gam',
                      eem_ind=FALSE,k=1,Rprobs=Rprobs)

  # evm results
  drcmd_evm_ests[ss] <- drcmd_res_evm$results$estimates$psi_hat_ate
  drcmd_direct_ests[ss] <- drcmd_res_evm$results$estimates$psi_hat_ate_direct
  temp_int <- compute_CI(drcmd_res_evm$results$estimates$psi_hat_ate,
                         drcmd_res_evm$results$ses$psi_hat_ate)
  drcmd_evm_cov[ss] <- (temp_int[1] < 2.5) & (temp_int[2] > 2.5)

  # mid results
  drcmd_mid_ests[ss] <- drcmd_res_mid$results$estimates$psi_hat_ate
  temp_int <- compute_CI(drcmd_res_mid$results$estimates$psi_hat_ate,
                         drcmd_res_mid$results$ses$psi_hat_ate)
  drcmd_mid_cov[ss] <- (temp_int[1] < 2.5) & (temp_int[2] > 2.5)
}



