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
