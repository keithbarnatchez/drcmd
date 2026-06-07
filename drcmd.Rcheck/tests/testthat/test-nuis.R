test_that("est_m_a returns predictions in [0, 1] for binary Y", {
  n <- 100
  X <- data.frame(x1 = rnorm(n))
  A <- rbinom(n, 1, plogis(X$x1))
  R <- rep(1, n)
  Y <- rbinom(n, 1, plogis(X$x1 + A))
  kappa_hat <- rep(1, n)
  idx <- 1:n
  m_learners <- c("SL.glm")

  out <- est_m_a(idx, Y, A, X, R, kappa_hat, m_learners)

  expect_true(all(out$m_1_hat >= 0 & out$m_1_hat <= 1))
  expect_true(all(out$m_0_hat >= 0 & out$m_0_hat <= 1))
})

test_that("est_g returns predictions in (0,1)", {
  n <- 200
  X <- data.frame(x1 = rnorm(n))
  A <- rbinom(n, 1, plogis(X$x1))
  R <- rep(1, n)
  kappa_hat <- rep(1, n)
  idx <- 1:n
  g_learners <- c("SL.glm")

  out <- est_g(idx, A, X, R, kappa_hat, g_learners)

  expect_true(all(out > 0 & out < 1))
})

test_that("est_kappa returns predictions in (0,1)", {
  n <- 200
  Z <- data.frame(z1 = rnorm(n))
  R <- rbinom(n, 1, plogis(Z$z1))
  idx <- 1:n
  r_learners <- c("SL.glm")

  out <- est_kappa(idx, Z, R, r_learners)

  expect_true(all(out > 0 & out < 1))
})

test_that("get_nuisance_ests returns all components", {
  n <- 200
  X <- data.frame(x1 = rnorm(n))
  A <- rbinom(n, 1, plogis(X$x1))
  Y <- A + X$x1 + rnorm(n)
  R <- rbinom(n, 1, plogis(X$x1))
  Y[R==0] <- 0
  Z <- X
  idx <- 1:n

  out <- get_nuisance_ests(idx, Y, A, X, Z, R,
                           'SL.glm', 'SL.glm', 'SL.glm',
                           Rprobs=NA, cutoff=0.025)

  expect_true("kappa_hat" %in% names(out))
  expect_true("m_a_hat" %in% names(out))
  expect_true("g_hat" %in% names(out))
  expect_equal(length(out$kappa_hat), n)
  expect_equal(length(out$g_hat), n)
})

# --- est_varphi tests ---

test_that("est_varphi returns correct structure for continuous Y", {
  set.seed(501)
  n <- 200
  Z <- data.frame(z1 = rnorm(n))
  R <- rbinom(n, 1, plogis(Z$z1))
  phi_1_hat <- rnorm(n)
  phi_0_hat <- rnorm(n)
  Y <- rnorm(n)
  idx <- 1:n

  out <- est_varphi(idx, R, Z, phi_1_hat, phi_0_hat,
                    po_learners = "SL.glm", Y = Y)

  expect_true("varphi_1_hat" %in% names(out))
  expect_true("varphi_0_hat" %in% names(out))
  expect_true("varphi_diff_hat" %in% names(out))
  expect_equal(length(out$varphi_1_hat), n)
  expect_equal(length(out$varphi_0_hat), n)
  expect_equal(length(out$varphi_diff_hat), n)
})

test_that("est_varphi handles binary Y (rescaling path)", {
  set.seed(502)
  n <- 200
  Z <- data.frame(z1 = rnorm(n))
  R <- rbinom(n, 1, plogis(Z$z1))
  phi_1_hat <- runif(n, 0.1, 0.9)
  phi_0_hat <- runif(n, 0.1, 0.9)
  Y <- rbinom(n, 1, 0.5)
  idx <- 1:n

  out <- est_varphi(idx, R, Z, phi_1_hat, phi_0_hat,
                    po_learners = "SL.glm", Y = Y)

  expect_equal(length(out$varphi_1_hat), n)
  expect_equal(length(out$varphi_0_hat), n)
})

# --- est_varphi_eem tests ---

test_that("est_varphi_eem returns correct structure", {
  set.seed(503)
  n <- 200
  Z <- data.frame(z1 = rnorm(n))
  R <- rbinom(n, 1, plogis(Z$z1))
  # Ensure kappa_hat is bounded away from 0 so R/kappa_hat - 1 doesn't blow up
  kappa_hat <- plogis(Z$z1)
  kappa_hat <- pmax(kappa_hat, 0.1)
  phi_1_hat <- rnorm(n)
  phi_0_hat <- rnorm(n)
  Y <- rnorm(n)
  idx <- 1:n

  out <- est_varphi_eem(idx, R, Z, phi_1_hat, phi_0_hat,
                        kappa_hat, po_learners = "SL.glm", Y = Y)

  expect_true("varphi_1_hat" %in% names(out))
  expect_true("varphi_0_hat" %in% names(out))
  expect_true("varphi_diff_hat" %in% names(out))
  expect_equal(length(out$varphi_1_hat), n)
})

# --- est_varphi_main tests ---

test_that("est_varphi_main returns zeros when R is all 1", {
  n <- 100
  R <- rep(1, n)
  Z <- data.frame(z1 = rnorm(n))
  kappa_hat <- rep(1, n)
  phi_1_hat <- rnorm(n)
  phi_0_hat <- rnorm(n)

  out <- est_varphi_main(1:n, R, Z, phi_1_hat, phi_0_hat,
                         kappa_hat, eem_ind = FALSE,
                         po_learners = "SL.glm", Y = rnorm(n))

  expect_true(all(out$varphi_1_hat == 0))
  expect_true(all(out$varphi_0_hat == 0))
  expect_true(all(out$varphi_diff_hat == 0))
})

test_that("est_varphi_main returns ATT/ATC varphi when requested (no missing data)", {
  n <- 100
  R <- rep(1, n)
  Z <- data.frame(z1 = rnorm(n))
  kappa_hat <- rep(1, n)
  phi_1_hat <- rnorm(n)
  phi_0_hat <- rnorm(n)
  phi_att_hat <- rnorm(n)
  phi_atc_hat <- rnorm(n)

  out <- est_varphi_main(1:n, R, Z, phi_1_hat, phi_0_hat,
                         kappa_hat, eem_ind = FALSE,
                         po_learners = "SL.glm", Y = rnorm(n),
                         att = TRUE, atc = TRUE,
                         phi_att_hat = phi_att_hat,
                         phi_atc_hat = phi_atc_hat)

  expect_true("varphi_att_hat" %in% names(out))
  expect_true("varphi_atc_hat" %in% names(out))
  expect_true(all(out$varphi_att_hat == 0))
  expect_true(all(out$varphi_atc_hat == 0))
})

test_that("est_varphi_main delegates to est_varphi with missing data", {
  set.seed(504)
  n <- 200
  Z <- data.frame(z1 = rnorm(n))
  R <- rbinom(n, 1, plogis(Z$z1))
  kappa_hat <- plogis(Z$z1)
  kappa_hat <- pmax(kappa_hat, 0.1)
  phi_1_hat <- rnorm(n)
  phi_0_hat <- rnorm(n)
  Y <- rnorm(n)

  out <- est_varphi_main(1:n, R, Z, phi_1_hat, phi_0_hat,
                         kappa_hat, eem_ind = FALSE,
                         po_learners = "SL.glm", Y = Y)

  expect_true("varphi_1_hat" %in% names(out))
  expect_equal(length(out$varphi_1_hat), n)
})

# --- SL.hal9001 / predict.SL.hal9001 tests ---

test_that("SL.hal9001 wrapper fits and predicts", {
  skip_if_not_installed("hal9001")
  set.seed(505)
  n <- 100
  X <- data.frame(x1 = rnorm(n))
  Y <- X$x1 + rnorm(n, sd = 0.5)

  fit <- SL.hal9001(Y, X, newX = X, family = gaussian(),
                    obsWeights = rep(1, n))

  expect_true("pred" %in% names(fit))
  expect_true("fit" %in% names(fit))
  expect_equal(length(fit$pred), n)
  expect_s3_class(fit$fit, "SL.hal9001")
})

test_that("predict.SL.hal9001 returns correct length", {
  skip_if_not_installed("hal9001")
  set.seed(506)
  n <- 100
  X <- data.frame(x1 = rnorm(n))
  Y <- X$x1 + rnorm(n, sd = 0.5)

  fit <- SL.hal9001(Y, X, newX = X, family = gaussian(),
                    obsWeights = rep(1, n))

  newX <- data.frame(x1 = rnorm(20))
  preds <- predict(fit$fit, newdata = newX)
  expect_equal(length(preds), 20)
})
