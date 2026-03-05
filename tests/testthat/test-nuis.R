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
