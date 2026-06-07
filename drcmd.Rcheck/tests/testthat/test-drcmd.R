context('Main functions of drcmd package')
set.seed(84123)


n <- 3000
X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
Y <- A + X + rnorm(n)/10 # small var for testing
Ystar <- Y + rnorm(n)/2
R <- rbinom(n,1,0.5*plogis(X)) # error-prone outcome measurements

# Make Y NA if R==0
Y[R==0] <- NA
X <- as.data.frame(X)

test_that("drcmd works with default parameters", {

  results <- drcmd(Y,A,X,
                   default_learners = 'SL.glm')

  expect_s3_class(results,"drcmd")

})

test_that('drcmd throws error for non-recognized SL libraries', {

  expect_error(drcmd(Y,A,X,
                     default_learners = 'SL.notanactuallibrary'))
})

test_that("one-step and tml give similar results", {

  results_onestep <- drcmd(Y,A,X,
                             default_learners = 'SL.glm')

  results_tml <- drcmd(Y,A,X,
                        default_learners = 'SL.glm',
                        tml = TRUE)

  expect_equal(results_onestep$results$estimates$psi_hat_ate,
               results_tml$results$estimates$psi_hat_ate,
               tolerance = 0.025)

})

test_that("drcmd works with binary outcome and reports RR/OR", {

  set.seed(7712)
  n <- 3000
  X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
  Y <- rbinom(n,1,plogis(X + A))
  R <- rbinom(n,1,0.5*plogis(X))
  Y[R==0] <- NA
  X <- as.data.frame(X)

  results <- drcmd(Y,A,X,
                   default_learners = 'SL.glm')

  expect_s3_class(results,"drcmd")
  expect_false(is.na(results$results$estimates$psi_hat_rr))
  expect_false(is.na(results$results$estimates$psi_hat_or))
  expect_false(is.na(results$results$ses$psi_hat_rr))
  expect_false(is.na(results$results$ses$psi_hat_or))

})

test_that("drcmd works with no missing data", {

  set.seed(9123)
  n <- 3000
  X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
  Y <- A + X + rnorm(n)/10
  X <- as.data.frame(X)

  results <- drcmd(Y,A,X,
                   default_learners = 'SL.glm')

  expect_s3_class(results,"drcmd")
  expect_true(all(results$R == 1))

})

test_that("drcmd works with cross-fitting (k>1)", {

  set.seed(5531)
  n <- 3000
  X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
  Y <- A + X + rnorm(n)/10
  R <- rbinom(n,1,0.5*plogis(X))
  Y[R==0] <- NA
  X <- as.data.frame(X)

  results_k1 <- drcmd(Y,A,X,
                       default_learners = 'SL.glm', k = 1)
  results_k2 <- drcmd(Y,A,X,
                       default_learners = 'SL.glm', k = 2)

  expect_s3_class(results_k2,"drcmd")
  expect_equal(results_k1$results$estimates$psi_hat_ate,
               results_k2$results$estimates$psi_hat_ate,
               tolerance = 0.1)

})

test_that("drcmd works with user-supplied Rprobs", {

  set.seed(4410)
  n <- 3000
  X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
  Y <- A + X + rnorm(n)/10
  rprobs <- plogis(X)
  R <- rbinom(n,1,rprobs)
  Y[R==0] <- NA
  X <- as.data.frame(X)

  results <- drcmd(Y,A,X,
                   default_learners = 'SL.glm',
                   Rprobs = rprobs)

  expect_s3_class(results,"drcmd")

})

test_that("drcmd works with multiple covariates", {

  set.seed(3319)
  n <- 3000
  X1 <- rnorm(n) ; X2 <- rnorm(n)
  A <- rbinom(n,1,plogis(X1 + X2))
  Y <- A + X1 - X2 + rnorm(n)/10
  R <- rbinom(n,1,plogis(X1))
  Y[R==0] <- NA
  X <- data.frame(X1=X1,X2=X2)

  results <- drcmd(Y,A,X,
                   default_learners = 'SL.glm')

  expect_s3_class(results,"drcmd")

})

test_that("ATE estimate is close to truth in well-specified case", {

  set.seed(1234)
  n <- 5000
  X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
  Y <- A + X + rnorm(n)/10 # true ATE = 1
  R <- rbinom(n,1,plogis(X))
  Y[R==0] <- NA
  X <- as.data.frame(X)

  results <- drcmd(Y,A,X,
                   default_learners = 'SL.glm')

  expect_equal(results$results$estimates$psi_hat_ate, 1,
               tolerance = 0.15)

})

test_that("returned object has expected structure", {

  set.seed(84123)
  n <- 3000
  X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
  Y <- A + X + rnorm(n)/10
  R <- rbinom(n,1,0.5*plogis(X))
  Y[R==0] <- NA
  X <- as.data.frame(X)

  results <- drcmd(Y,A,X,
                   default_learners = 'SL.glm')

  expect_true("results" %in% names(results))
  expect_true("params" %in% names(results))
  expect_true("Z" %in% names(results))
  expect_true("U" %in% names(results))
  expect_true("R" %in% names(results))

  expect_true("estimates" %in% names(results$results))
  expect_true("ses" %in% names(results$results))
  expect_true("nuis" %in% names(results$results))

  expect_true(all(c("psi_1_hat","psi_0_hat","psi_hat_ate") %in%
                    colnames(results$results$estimates)))

})

# --- ATT / ATC tests ---

test_that("ATT/ATC are NA by default", {

  set.seed(1111)
  n <- 2000
  X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
  Y <- A + X + rnorm(n)/10
  R <- rbinom(n,1,plogis(X))
  Y[R==0] <- NA
  X <- as.data.frame(X)

  results <- drcmd(Y,A,X, default_learners='SL.glm')

  expect_true(is.na(results$results$estimates$psi_hat_att))
  expect_true(is.na(results$results$estimates$psi_hat_atc))
  expect_true(is.na(results$results$ses$psi_hat_att))
  expect_true(is.na(results$results$ses$psi_hat_atc))

})

test_that("ATT estimate is close to truth", {

  set.seed(2222)
  n <- 5000
  X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
  Y <- A + X + rnorm(n)/10 # true ATT = 1 (constant treatment effect)
  R <- rbinom(n,1,plogis(X))
  Y[R==0] <- NA
  X <- as.data.frame(X)

  results <- drcmd(Y,A,X, default_learners='SL.glm', att=TRUE)

  expect_false(is.na(results$results$estimates$psi_hat_att))
  expect_equal(results$results$estimates$psi_hat_att, 1, tolerance=0.15)

})

test_that("ATC estimate is close to truth", {

  set.seed(3333)
  n <- 5000
  X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
  Y <- A + X + rnorm(n)/10 # true ATC = 1 (constant treatment effect)
  R <- rbinom(n,1,plogis(X))
  Y[R==0] <- NA
  X <- as.data.frame(X)

  results <- drcmd(Y,A,X, default_learners='SL.glm', atc=TRUE)

  expect_false(is.na(results$results$estimates$psi_hat_atc))
  expect_equal(results$results$estimates$psi_hat_atc, 1, tolerance=0.15)

})

test_that("ATT/ATC work with cross-fitting (k>1)", {

  set.seed(4444)
  n <- 3000
  X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
  Y <- A + X + rnorm(n)/10
  R <- rbinom(n,1,plogis(X))
  Y[R==0] <- NA
  X <- as.data.frame(X)

  results <- drcmd(Y,A,X, default_learners='SL.glm', k=2,
                   att=TRUE, atc=TRUE)

  expect_false(is.na(results$results$estimates$psi_hat_att))
  expect_false(is.na(results$results$estimates$psi_hat_atc))
  expect_false(is.na(results$results$ses$psi_hat_att))
  expect_false(is.na(results$results$ses$psi_hat_atc))

})

test_that("ATT/ATC are NA under TML path", {

  set.seed(5555)
  n <- 3000
  X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
  Y <- rbinom(n,1,plogis(X + A))
  R <- rbinom(n,1,plogis(X))
  Y[R==0] <- NA
  X <- as.data.frame(X)

  results <- drcmd(Y,A,X, default_learners='SL.glm',
                   tml=TRUE, att=TRUE, atc=TRUE)

  expect_true(is.na(results$results$estimates$psi_hat_att))
  expect_true(is.na(results$results$estimates$psi_hat_atc))

})

test_that("ATT/ATC work with no missing data", {

  set.seed(6666)
  n <- 3000
  X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
  Y <- A + X + rnorm(n)/10
  X <- as.data.frame(X)

  results <- drcmd(Y,A,X, default_learners='SL.glm',
                   att=TRUE, atc=TRUE)

  expect_false(is.na(results$results$estimates$psi_hat_att))
  expect_false(is.na(results$results$estimates$psi_hat_atc))
  expect_equal(results$results$estimates$psi_hat_att, 1, tolerance=0.15)
  expect_equal(results$results$estimates$psi_hat_atc, 1, tolerance=0.15)

})

# --- get_phi_hat tests ---

test_that("get_phi_hat returns correct structure and length", {

  set.seed(7001)
  n <- 200
  X <- data.frame(x1 = rnorm(n))
  A <- rbinom(n, 1, 0.5)
  Y <- A + X$x1 + rnorm(n)
  R <- rep(1, n)
  Z <- X
  g_hat <- rep(0.5, n)
  m_a_hat <- list(m_1_hat = 1 + X$x1, m_0_hat = X$x1)
  kappa_hat <- rep(1, n)

  phi <- get_phi_hat(Y, A, X, R, Z, g_hat, m_a_hat, kappa_hat)

  expect_true("phi_1_hat" %in% names(phi))
  expect_true("phi_0_hat" %in% names(phi))
  expect_true("plugin1" %in% names(phi))
  expect_true("plugin0" %in% names(phi))
  expect_equal(length(phi$phi_1_hat), n)
  expect_equal(length(phi$phi_0_hat), n)

})

test_that("get_phi_hat uses IPW plugin when X not subset of Z", {

  set.seed(7002)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Z <- data.frame(x1 = X$x1) # Z only has x1, not x2
  A <- rbinom(n, 1, 0.5)
  Y <- rnorm(n)
  R <- rbinom(n, 1, 0.7)
  g_hat <- rep(0.5, n)
  m_a_hat <- list(m_1_hat = rep(1, n), m_0_hat = rep(0, n))
  kappa_hat <- rep(0.7, n)

  phi <- get_phi_hat(Y, A, X, R, Z, g_hat, m_a_hat, kappa_hat)

  # When X not in Z, plugin uses IPW: mean(R/kappa * m_hat)
  expect_equal(phi$plugin1, mean(R / kappa_hat * m_a_hat$m_1_hat))

})

# --- tml_updates tests ---

test_that("tml_updates returns updated nuisance estimates", {

  set.seed(7003)
  n <- 500
  X <- data.frame(x1 = rnorm(n))
  A <- rbinom(n, 1, plogis(X$x1))
  Y <- rbinom(n, 1, plogis(X$x1 + A))
  R <- rbinom(n, 1, plogis(X$x1))
  Z <- X

  # Fit nuisance models
  nuis <- get_nuisance_ests(1:n, Y, A, X, Z, R,
                            'SL.glm', 'SL.glm', 'SL.glm',
                            Rprobs = NA, cutoff = 0.025)

  phi_hat <- get_phi_hat(Y, A, X, R, Z,
                         nuis$g_hat, nuis$m_a_hat, nuis$kappa_hat)

  varphi_hat <- est_varphi_main(1:n, R, Z,
                                phi_hat$phi_1_hat, phi_hat$phi_0_hat,
                                nuis$kappa_hat, eem_ind = FALSE,
                                po_learners = "SL.glm", Y = Y)

  updated <- tml_updates(1:n, Y, A, X, R, Z,
                         nuis$m_a_hat$m_1_hat, nuis$m_a_hat$m_0_hat,
                         nuis$g_hat, nuis$kappa_hat,
                         phi_hat$phi_1_hat, phi_hat$phi_0_hat,
                         varphi_hat)

  expect_true("m_1_hat_star" %in% names(updated))
  expect_true("m_0_hat_star" %in% names(updated))
  expect_true("plugin_ate_star" %in% names(updated))
  expect_true("kappa_hat_ate_star" %in% names(updated))
  expect_equal(length(updated$m_1_hat_star), n)
  # Updated m predictions should be in (0,1) since Y is binary
  expect_true(all(updated$m_1_hat_star >= 0 & updated$m_1_hat_star <= 1))
  expect_true(all(updated$m_0_hat_star >= 0 & updated$m_0_hat_star <= 1))

})
