# test_funcs.R
# space for testing functions in the drcmd package
#-------------------------------------------------------------------------------
source('utils.R')
source('drcmd.R')
#-------------------------------------------------------------------------------
# Params for functions

hal_ind <- FALSE
eem_ind <- FALSE
sl.lib <- 'SL.glm'
#-------------------------------------------------------------------------------
# Make a couple functions for simming simple data structure
n <- 1e3
X <- rnorm(n) ; A <- rbinom(n,1,plogis(X)) ; Y <- rnorm(n) + A + X
Ystar <- Y + rnorm(n)/2 ; R <- rbinom(n,1,plogis(X))

# Make Y NA is R==0
Y[R==0] <- NA

df <- data.frame(Y=Y,A=A,X=X,Ystar=Ystar,R=R)

X <- 'X' ; A <- 'A' ; Y <- 'Y' ; R <- 'R' ; W <- 'Ystar'
#-------------------------------------------------------------------------------
# Test sub functions

V <- find_missing_pattern(df,Y,A,X,W,R)
Z <- V$Z ; R <- V$Rstr ; df[,R] <- V$Rvals

kappa_hat <- est_kappa(df,Z,R,hal_ind,sl.lib)
m_hat <- est_m_a(df,Y,A,X,R,kappa_hat,hal_ind,sl.lib)
g_hat <- est_g(df,A,X,R,kappa_hat,hal_ind,sl.lib)

# put together those jawns
m_1_hat <- m_hat$m_1_hat
m_0_hat <- m_hat$m_0_hat
# g_hat <- g_hat$g_hats

# Check the phi maker
phi_hat  <- get_phi_hat(df,Y,A,X,Z,R,g_hat,m_hat,kappa_hat,hal_ind,sl.lib)
phi_1_hat <- phis$phi_1_hat ; phi_0_hat <- phis$phi_0_hat

varphi_hat <- est_varphi(df,R,Z,phi_1_hat,phi_0_hat,hal_ind,sl.lib)

est_psi(df, R, Z, kappa_hat, phi_hat,varphi_hat)
#-------------------------------------------------------------------------------
# Test outer functions
temp <- drcmd_est_fold(df,splits,Y,A,X,Z,R,hal_ind,sl.lib,eem_ind,Rprobs)

splits <- create_folds(df,10)
temp <- drcmd_est(df,Y,A,X,Z,R,hal_ind,sl.lib,eem_ind,Rprobs,k=3)
