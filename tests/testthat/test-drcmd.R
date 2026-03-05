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
