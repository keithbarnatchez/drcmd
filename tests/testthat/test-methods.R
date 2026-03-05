# Simulate data to pass into drcmd
n <- 3000
X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
Y <- A + X + rnorm(n)/10 # small var for testing
Ystar <- Y + rnorm(n)/2
R <- rbinom(n,1,0.5*plogis(X)) # error-prone outcome measurements

# Make Y NA if R==0
Y[R==0] <- NA
X <- as.data.frame(X)

results <- drcmd(Y,A,X,
                 default_learners = 'SL.glm')

test_that('error thrown if non-standard plot type provided', {
  expect_error(plot(results,type='not a valid type'))
})

test_that('error thrown if wrong arg thrown at summary', {
  expect_error(summary(results,detail='this should be a logical'))
})

test_that("print.drcmd runs without error", {
  expect_output(print(results), "drcmd results")
  expect_output(print(results), "ATE estimate")
})

test_that("summary.drcmd runs without error", {
  expect_output(summary(results), "Summary of drcmd results")
  expect_output(summary(results), "ATE")
})

test_that("summary.drcmd with detail=TRUE prints learner info", {
  expect_output(summary(results, detail=TRUE), "nuisance learners")
})
