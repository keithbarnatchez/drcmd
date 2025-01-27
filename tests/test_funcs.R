# test_funcs.R
# space for testing functions in the drcmd package
#-------------------------------------------------------------------------------
rm(list=ls())
source('../R/utils.R')
source('../R/drcmd.R')
source('../R/nuis.R')
#-------------------------------------------------------------------------------
# Params for functions

hal_ind <- FALSE
eem_ind <- FALSE
sl.lib <- 'SL.glm'
#-------------------------------------------------------------------------------
# Make a couple functions for simming simple data structure
n <- 1e3
X <- rnorm(n) ; A <- rbinom(n,1,plogis(X)) ; Y <- rnorm(n) + A + X
Ystar <- Y + rnorm(n)/2 ; R <- rbinom(n,1,plogis(X)) ; X <- as.data.frame(X)


# Make Y NA is R==0
Y[R==0] <- NA

df <- data.frame(Y=Y,A=A,X=X,Ystar=Ystar,R=R)

# X = data.frame(cbind(X,X))
W = X[,0]

drcmd_res <- drcmd(Y,A,X, hal_ind=FALSE,sl.lib=sl.lib)
#-------------------------------------------------------------------------------
# Test sub functions

V <- find_missing_pattern(Y,A,X,W)
Z <- V$Z ; R <- V$R ; X <- V$X ; Y <- V$Y ; A <- V$A
Rprobs <- NA

# nuisance functions
idx <- 1:n
kappa_hat <- est_kappa(idx,Z,R,hal_ind,sl.lib)
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

