# Package index

## Estimation

Main entry point for doubly-robust estimation with missing data.

- [`drcmd()`](https://kbarnatchez.com/drcmd/reference/drcmd.md) :
  Doubly-robust causal inference with missing data

## Nuisance estimation

Helpers for specifying and fitting nuisance functions.

- [`get_sl_libraries()`](https://kbarnatchez.com/drcmd/reference/get_sl_libraries.md)
  : List SuperLearner libraries
- [`SL.hal9001()`](https://kbarnatchez.com/drcmd/reference/SL.hal9001.md)
  : SuperLearner wrapper for the highly-adaptive lasso

## Methods

S3 methods for drcmd objects.

- [`plot(`*`<drcmd>`*`)`](https://kbarnatchez.com/drcmd/reference/plot.drcmd.md)
  : Plot results from drcmd object
- [`print(`*`<drcmd>`*`)`](https://kbarnatchez.com/drcmd/reference/print.drcmd.md)
  : Print drcmd object
- [`summary(`*`<drcmd>`*`)`](https://kbarnatchez.com/drcmd/reference/summary.drcmd.md)
  : Summarize results from drcmd
- [`predict(`*`<SL.hal9001>`*`)`](https://kbarnatchez.com/drcmd/reference/predict.SL.hal9001.md)
  : Prediction wrapper for the highly-adaptive lasso
