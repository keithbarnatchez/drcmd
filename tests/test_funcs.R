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
default_learners <- 'SL.glm'
#-------------------------------------------------------------------------------
# Make a couple functions for simming simple data structure
n <- 1e3
X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
Y <- rnorm(n) + A + X + X^2 + A*X + sin(X)
Ystar <- Y + rnorm(n)/2 ; R <- rbinom(n,1,plogis(X)) ; X <- as.data.frame(X)

# Make Y NA if R==0
Y[R==0] <- NA

df <- data.frame(Y=Y,A=A,X=X,Ystar=Ystar,R=R)

# X = data.frame(cbind(X,X))
W = X[,0]

drcmd_res <- drcmd(Y,A,X, default_learners= c('SL.glm','SL.glm.interaction','SL.earth','SL.gam'),
                   eem_ind=FALSE,k=1)
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

df <- data.frame(Y=Y,A=A,X=X,Ystar=Ystar,R=R)

W <- data.frame(Ystar)

drcmd_res <- drcmd(Y,A,X,
                   W=data.frame(Ystar),
                   default_learners= c('SL.glm','SL.glm.interaction','SL.earth','SL.gam'),
                   eem_ind=FALSE,k=1)
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

