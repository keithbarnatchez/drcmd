# Simulation tests can legitimately trigger truncation or select zero
# augmentation. Check each such diagnostic explicitly; never muffle other
# warnings (in particular convergence, rank deficiency, or invalid fits).
expect_fit_diagnostics <- function(expr, truncation = FALSE, augmentation = FALSE) {
  withCallingHandlers(expr, warning = function(w) {
    if (augmentation && inherits(w, "drcmd_zero_augmentation")) {
      expect_match(conditionMessage(w), "predictions are zero", fixed = TRUE)
      invokeRestart("muffleWarning")
    }
    if (truncation && grepl(
      "^(Propensity scores|Complete case probabilities) outside of [0-9.]+ and [0-9.]+\\. Truncating to cutoffs$",
      conditionMessage(w))) {
      expect_match(conditionMessage(w), "Truncating to cutoffs", fixed = TRUE)
      invokeRestart("muffleWarning")
    }
  })
}

drcmd_test_fit <- function(...) {
  expect_fit_diagnostics(drcmd(...), truncation = TRUE, augmentation = TRUE)
}
