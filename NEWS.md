# drcmd 0.1.0

* Reject degenerate treatment and complete-case probability fits before
  truncation, with a message identifying the failed nuisance regression.
* Skip unnecessary TML missingness updates for fully observed data and genuine
  zero-clever-covariate updates. Stop on nonconverged, unidentified, or nonfinite
  targeting coefficients instead of silently treating the update as zero.
* Use quasibinomial outcome targeting to avoid spurious binomial-count warnings
  for fractional outcomes and inverse-probability weights.
* Report zero augmentation ensembles with an explicit warning rather than
  silently suppressing the diagnostic.
* Remove the ignored `reduce_basis` argument from the optional HAL wrapper,
  retaining the underlying learner's existing smoothness default.
* Add regression tests for degenerate fits and TML updates, and explicitly test
  expected fitting diagnostics without suppressing unrelated warnings.
