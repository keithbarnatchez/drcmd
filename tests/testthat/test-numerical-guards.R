test_that("probability fits reject zero ensembles without hiding other warnings", {
  for (msg in c("All algorithms have zero weight",
                "All metalearner coefficients are zero, predictions will all be equal to 0")) {
    expect_error(with_probability_fit_checks(warning(msg), "Treatment regression"),
                 "Treatment regression produced an all-zero")
  }
  expect_silent(with_probability_fit_checks({
    warning("non-integer #successes in a binomial glm!")
    1
  }, "Treatment regression"))
  expect_warning(with_probability_fit_checks(
    warning("glm.fit: algorithm did not converge"), "Treatment regression"),
    "algorithm did not converge")
})

test_that("probability predictions are checked before truncation", {
  for (x in list(numeric(), c(NA, 0.5), c(Inf, 0.5), c(-0.1, 0.5), c(0.5, 1.1))) {
    expect_error(check_probability_predictions(x, "Complete-case regression"),
                 "invalid probability predictions")
  }
  for (x in list(rep(0, 10), rep(1, 10))) {
    expect_error(check_probability_predictions(x, "Treatment regression"),
                 "only boundary probabilities")
  }
  # Ordinary constant probabilities and individual boundary predictions are
  # allowed; the latter remain subject to the user's truncation choice.
  expect_silent(check_probability_predictions(rep(0.5, 10), "Treatment regression"))
  expect_silent(check_probability_predictions(c(0, 0.5, 1), "Treatment regression"))
})

test_that("zero augmentation is allowed with an explicit diagnostic", {
  set.seed(701)
  X <- data.frame(x = seq_len(40))
  expect_warning(fit <- fit_augmentation_sl(
    Y = rep(c(-1, 1), each = 20), X = X, family = gaussian(), SL.library = "SL.mean",
    cvControl = list(V = 2, validRows = list(1:20, 21:40))),
    class = "drcmd_zero_augmentation")
  expect_true(all(fit$coef == 0))
  expect_silent(pred <- predict_augmentation_sl(fit, X))
  expect_equal(pred, rep(0, 40))

  expect_silent(nonzero <- fit_augmentation_sl(
    Y = 2 + X$x, X = X, family = gaussian(), SL.library = "SL.glm",
    cvControl = list(V = 2)))
  expect_equal(as.vector(predict_augmentation_sl(nonzero, X)), 2 + X$x,
               tolerance = 1e-8)
})

test_that("failed TML updates cannot silently become no-ops", {
  valid <- list(coefficients = c(epsilon = 0.2), converged = TRUE, rank = 1L)
  expect_equal(targeting_coefficient(valid, "ate", "outcome"), 0.2)
  for (bad in list(
    modifyList(valid, list(converged = FALSE)),
    modifyList(valid, list(coefficients = NA_real_)),
    modifyList(valid, list(coefficients = Inf)),
    modifyList(valid, list(coefficients = numeric())),
    modifyList(valid, list(rank = 0L)))) {
    expect_error(targeting_coefficient(bad, "ate", "outcome"),
                 "TML ate outcome update failed")
  }
})

test_that("TML handles full observation and zero clever covariates exactly", {
  n <- 80
  A <- rep(c(0, 1), n / 2)
  Y <- plogis(seq(-1, 1, length.out = n) + 0.2 * A)
  args <- list(idx = seq_len(n), Y = Y, A = A, R = rep(1, n),
               g_hat = rep(0.5, n), kappa_hat = rep(1, n),
               m_1_hat = rep(0.6, n), m_0_hat = rep(0.4, n),
               varphi_1 = rep(0, n), varphi_0 = rep(0, n), varphi_diff = rep(0, n))
  for (target in c("ate", "psi_1", "psi_0")) {
    expect_silent(out <- do.call(.tml_target_one, c(list(target = target), args)))
    expect_identical(out$kappa_hat_star, rep(1, n))
    expect_true(all(is.finite(out$m_1_hat_star)))
    expect_true(all(is.finite(out$m_0_hat_star)))
  }

  # No treated observations in this evaluation fold: no outcome update for psi_1.
  args$idx <- which(A == 0)
  expect_silent(out <- do.call(.tml_target_one, c(list(target = "psi_1"), args)))
  expect_equal(out$m_1_hat_star, args$m_1_hat)
  expect_equal(out$m_0_hat_star, args$m_0_hat)

  # Zero augmentation also means no missingness fluctuation when R varies.
  args$idx <- seq_len(n)
  args$R <- rep(c(1, 1, 0, 1), length.out = n)
  args$kappa_hat <- rep(0.7, n)
  expect_silent(out <- do.call(.tml_target_one, c(list(target = "ate"), args)))
  expect_identical(out$kappa_hat_star, args$kappa_hat)

  # Changing binomial to quasibinomial must preserve the targeting coefficient.
  H <- A / args$g_hat - (1 - A) / (1 - args$g_hat)
  offset <- qlogis(args$m_1_hat * A + args$m_0_hat * (1 - A))
  expect_warning(old <- glm(Y ~ -1 + H, offset = offset,
                            weights = args$R / args$kappa_hat, family = binomial()),
                 "non-integer #successes")
  expect_equal(out$m_1_hat_star,
               plogis(qlogis(args$m_1_hat) + as.double(coef(old)) / args$g_hat))
  expect_equal(out$m_0_hat_star,
               plogis(qlogis(args$m_0_hat) - as.double(coef(old)) / (1 - args$g_hat)))
})
