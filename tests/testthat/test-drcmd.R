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
