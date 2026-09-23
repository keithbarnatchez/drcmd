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

  results <- drcmd_test_fit(Y,A,X,
                   default_learners = 'SL.glm')

  expect_s3_class(results,"drcmd")

})

test_that('drcmd throws error for non-recognized SL libraries', {

  expect_error(drcmd_test_fit(Y,A,X,
                     default_learners = 'SL.notanactuallibrary'))
})

test_that("one-step and tml give similar results", {

  results_onestep <- drcmd_test_fit(Y,A,X,
                             default_learners = 'SL.glm')

  results_tml <- drcmd_test_fit(Y,A,X,
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

  results <- drcmd_test_fit(Y,A,X,
                   default_learners = 'SL.glm')

  expect_s3_class(results,"drcmd")
  expect_false(is.na(results$results$estimates$psi_hat_rr))
  expect_false(is.na(results$results$estimates$psi_hat_or))
  expect_false(is.na(results$results$ses$psi_hat_rr))
  expect_false(is.na(results$results$ses$psi_hat_or))

})

test_that("SL.ranger fits binary-outcome pseudo-outcome regressions", {

  skip_if_not_installed("ranger")

  set.seed(9106)
  n <- 240
  X <- data.frame(X = rnorm(n))
  A <- rbinom(n, 1, plogis(0.2 + 0.5 * X$X))
  Y <- rbinom(n, 1, plogis(-0.4 + 0.8 * A + 0.6 * X$X))
  R <- rbinom(n, 1, plogis(0.7 - 0.4 * X$X))
  Y[R == 0] <- NA

  for (k in c(1, 2)) {
    set.seed(9106)
    fit <- drcmd_test_fit(
      Y, A, X,
      m_learners = "SL.glm",
      g_learners = "SL.glm",
      r_learners = "SL.glm",
      po_learners = "SL.ranger",
      k = k,
      cv_folds = 2
    )

    expect_s3_class(fit, "drcmd")
    expect_true(all(is.finite(fit$results$nuis$varphi_1_hat)))
    expect_true(all(is.finite(fit$results$nuis$varphi_0_hat)))
  }
})

test_that("drcmd works with no missing data", {

  set.seed(9123)
  n <- 3000
  X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
  Y <- A + X + rnorm(n)/10
  X <- as.data.frame(X)

  results <- drcmd_test_fit(Y,A,X,
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

  results_k1 <- drcmd_test_fit(Y,A,X,
                       default_learners = 'SL.glm', k = 1)
  results_k2 <- drcmd_test_fit(Y,A,X,
                       default_learners = 'SL.glm', k = 2)

  expect_s3_class(results_k2,"drcmd")
  expect_equal(results_k1$results$estimates$psi_hat_ate,
               results_k2$results$estimates$psi_hat_ate,
               tolerance = 0.1)

})

test_that("cross-fitting identifies insufficient binary nuisance classes", {

  set.seed(5532)
  n <- 40
  X <- data.frame(X = rnorm(n))
  A <- c(1, rep(0, n - 1))
  Y <- A + X$X + rnorm(n)

  expect_warning(expect_error(
    drcmd_test_fit(Y, A, X, default_learners = "SL.glm", k = 2, cv_folds = 2),
    "Cross-fitting fold 1: Treatment regression requires at least two"
  ), "prediction from rank-deficient fit")
})

test_that("cross-fitted binary contrasts and standard errors are well formed", {

  set.seed(20260821)
  n <- 800
  X <- rnorm(n)
  A <- rbinom(n, 1, plogis(X))
  Y <- rbinom(n, 1, plogis(-0.5 + 0.8 * A + 0.5 * X))
  R <- rbinom(n, 1, plogis(0.6 + 0.3 * X))
  Y[R == 0] <- NA
  X <- data.frame(X = X)

  one_step <- drcmd_test_fit(Y, A, X, default_learners = "SL.glm",
                    k = 2, cv_folds = 2)
  set.seed(20260821)
  tml_fit <- drcmd_test_fit(Y, A, X, default_learners = "SL.glm",
                   k = 2, cv_folds = 2, tml = TRUE)

  for (fit in list(one_step, tml_fit)) {
    est <- fit$results$estimates
    se <- fit$results$ses

    expect_equal(est$psi_hat_rr, est$psi_1_hat / est$psi_0_hat)
    expect_equal(
      est$psi_hat_or,
      (est$psi_1_hat / (1 - est$psi_1_hat)) /
        (est$psi_0_hat / (1 - est$psi_0_hat))
    )
    expect_true(all(is.finite(unlist(se[c("psi_hat_rr", "psi_hat_or")]))))
    expect_lt(se$psi_hat_rr, 2)
    expect_lt(se$psi_hat_or, 2)
    expect_equal(nrow(fit$results$nuis), n)
  }

})

test_that("cross-fitted TML does not report binary contrasts for continuous outcomes", {

  set.seed(20260822)
  n <- 400
  X <- rnorm(n)
  A <- rbinom(n, 1, plogis(X))
  Y <- A + X + rnorm(n)
  X <- data.frame(X = X)

  fit <- drcmd_test_fit(Y, A, X, default_learners = "SL.glm",
               k = 2, cv_folds = 2, tml = TRUE)

  expect_true(is.na(fit$results$estimates$psi_hat_rr))
  expect_true(is.na(fit$results$estimates$psi_hat_or))
  expect_true(is.na(fit$results$ses$psi_hat_rr))
  expect_true(is.na(fit$results$ses$psi_hat_or))

})

test_that("TML fits bounded continuous outcomes with SL.ranger", {

  skip_if_not_installed("ranger")

  set.seed(20260906)
  n <- 240
  X <- data.frame(X = runif(n))
  A <- rbinom(n, 1, plogis(-0.2 + X$X))
  Y <- pmin(pmax(0.2 + 0.3 * A + 0.2 * X$X + rnorm(n, sd = 0.08), 0), 1)
  Y[1] <- 0
  Y[2] <- 1

  for (k in c(1, 2)) {
    set.seed(20260906)
    fit <- drcmd_test_fit(Y, A, X, default_learners = "SL.ranger",
                 k = k, cv_folds = 2, tml = TRUE)

    expect_s3_class(fit, "drcmd")
    expect_true(is.finite(fit$results$estimates$psi_hat_ate))
    expect_true(is.na(fit$results$estimates$psi_hat_rr))
  }
})

test_that("cross-fitted delta-method variances use point estimates", {

  fold_1 <- list(ics = data.frame(
    psi_1_ic = c(-0.2, 0.1),
    psi_0_ic = c(-0.1, 0.05),
    psi_ate_ic = c(-0.1, 0.05),
    psi_att_ic = NA_real_,
    psi_atc_ic = NA_real_
  ))
  fold_2 <- list(ics = data.frame(
    psi_1_ic = c(0.05, 0.05),
    psi_0_ic = c(0.02, 0.03),
    psi_ate_ic = c(0.03, 0.02),
    psi_att_ic = NA_real_,
    psi_atc_ic = NA_real_
  ))
  ests <- c(psi_1_hat = 0.6, psi_0_hat = 0.3)

  got <- est_ses_crossfit(list(fold_1, fold_2), ests, y_bin = TRUE)
  ic <- rbind(fold_1$ics, fold_2$ics)
  sig <- cov(ic[, c("psi_1_ic", "psi_0_ic")])
  expected_rr_var <- (sig[1, 1] / ests["psi_0_hat"]^2 -
                        2 * sig[1, 2] * ests["psi_1_hat"] /
                          ests["psi_0_hat"]^3 +
                        sig[2, 2] * ests["psi_1_hat"]^2 /
                          ests["psi_0_hat"]^4) / nrow(ic)

  expect_equal(got$psi_hat_rr, unname(expected_rr_var))

})

test_that("pseudo-outcome nuisance fitting does not use held-out outcomes", {

  set.seed(8201)
  n <- 400
  X <- data.frame(X = rnorm(n))
  A <- rbinom(n, 1, plogis(X$X))
  Y <- A + X$X + rnorm(n)
  R <- rbinom(n, 1, plogis(0.5 + X$X))
  Y[R == 0] <- 0
  Z <- X
  splits <- list(train = 1:300, test = 301:400)

  set.seed(8202)
  fit_1 <- expect_fit_diagnostics(drcmd_est_fold(
    splits, Y, A, X, Z, R,
    m_learners = "SL.glm", g_learners = "SL.glm",
    r_learners = "SL.glm", po_learners = "SL.glm",
    eem_ind = FALSE, tml = FALSE, Rprobs = NA,
    cutoff = 0.025, y_bin = FALSE, cv_folds = 2
  ), truncation = TRUE, augmentation = TRUE)

  Y_changed <- Y
  Y_changed[splits$test] <- Y_changed[splits$test] + 100
  set.seed(8202)
  fit_2 <- expect_fit_diagnostics(drcmd_est_fold(
    splits, Y_changed, A, X, Z, R,
    m_learners = "SL.glm", g_learners = "SL.glm",
    r_learners = "SL.glm", po_learners = "SL.glm",
    eem_ind = FALSE, tml = FALSE, Rprobs = NA,
    cutoff = 0.025, y_bin = FALSE, cv_folds = 2
  ), truncation = TRUE, augmentation = TRUE)

  expect_equal(fit_1$nuis$varphi_1_hat, fit_2$nuis$varphi_1_hat)
  expect_equal(fit_1$nuis$varphi_0_hat, fit_2$nuis$varphi_0_hat)

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

  results <- drcmd_test_fit(Y,A,X,
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

  results <- drcmd_test_fit(Y,A,X,
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

  results <- drcmd_test_fit(Y,A,X,
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

  results <- drcmd_test_fit(Y,A,X,
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

  results <- drcmd_test_fit(Y,A,X, default_learners='SL.glm')

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

  results <- drcmd_test_fit(Y,A,X, default_learners='SL.glm', att=TRUE)

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

  results <- drcmd_test_fit(Y,A,X, default_learners='SL.glm', atc=TRUE)

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

  results <- drcmd_test_fit(Y,A,X, default_learners='SL.glm', k=2,
                   att=TRUE, atc=TRUE)

  expect_false(is.na(results$results$estimates$psi_hat_att))
  expect_false(is.na(results$results$estimates$psi_hat_atc))
  expect_false(is.na(results$results$ses$psi_hat_att))
  expect_false(is.na(results$results$ses$psi_hat_atc))

})

test_that("ATT/ATC error under TML path", {

  set.seed(5555)
  n <- 3000
  X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
  Y <- rbinom(n,1,plogis(X + A))
  R <- rbinom(n,1,plogis(X))
  Y[R==0] <- NA
  X <- as.data.frame(X)

  # ATT/ATC are only supported under the one-step estimator
  expect_error(
    drcmd_test_fit(Y,A,X, default_learners='SL.glm',
          tml=TRUE, att=TRUE, atc=TRUE),
    "ATT/ATC estimation is only supported with one-step estimation"
  )

})

test_that("ATT/ATC work with no missing data", {

  set.seed(6666)
  n <- 3000
  X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
  Y <- A + X + rnorm(n)/10
  X <- as.data.frame(X)

  results <- drcmd_test_fit(Y,A,X, default_learners='SL.glm',
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

  varphi_hat <- expect_fit_diagnostics(est_varphi_main(1:n, R, Z,
                                phi_hat$phi_1_hat, phi_hat$phi_0_hat,
                                nuis$kappa_hat, eem_ind = FALSE,
                                po_learners = "SL.glm", Y = Y), augmentation = TRUE)

  updated <- tml_updates(1:n, Y, A, X, R, Z,
                         nuis$m_a_hat$m_1_hat, nuis$m_a_hat$m_0_hat,
                         nuis$g_hat, nuis$kappa_hat,
                         phi_hat$phi_1_hat, phi_hat$phi_0_hat,
                         varphi_hat)

  # tml_updates runs a separate fluctuation per estimand
  expect_true(all(c("ate", "psi_1", "psi_0") %in% names(updated)))
  expect_true("m_1_hat_star" %in% names(updated$ate))
  expect_true("m_0_hat_star" %in% names(updated$ate))
  expect_true("kappa_hat_star" %in% names(updated$ate))
  expect_equal(length(updated$ate$m_1_hat_star), n)
  # Updated m predictions should be in (0,1) since Y is binary
  expect_true(all(updated$ate$m_1_hat_star >= 0 & updated$ate$m_1_hat_star <= 1))
  expect_true(all(updated$ate$m_0_hat_star >= 0 & updated$ate$m_0_hat_star <= 1))

})
