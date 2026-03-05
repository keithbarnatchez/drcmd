context('Helper functions of drcmd package')

test_that('drcmd throws error if no complete cases', {
  n <- 3000
  X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
  Y <- A + X + rnorm(n)/10 # small var for testing
  Ystar <- Y + rnorm(n)/2

  # Make Y NA if R==0
  Y <- NA
  X <- as.data.frame(X)

  expect_error(drcmd(Y,A,X,
             default_learners = 'SL.glm'))
})

test_that('drcmd throws error if no columns are fully non-missing', {
  n <- 3000
  X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
  Y <- A + X + rnorm(n)/10 # small var for testing
  Ystar <- Y + rnorm(n)/2
  R <- rbinom(n,1,0.5*plogis(X)) # error-prone outcome measurements

  # Make Y NA if R==0
  Y[R==0] <- NA
  X[R==0] <- NA
  A[R==0] <- NA
  X <- as.data.frame(X)

  expect_error(drcmd(Y,A,X,
             default_learners = 'SL.glm'))
})

test_that('check_entry_errors catches issue upfront', {
  n <- 3000
  X <- rnorm(n) ; A <- rbinom(n,1,plogis(X))
  Y <- A + X + rnorm(n)/10 # small var for testing
  Ystar <- Y + rnorm(n)/2
  R <- rbinom(n,1,0.5*plogis(X)) # error-prone outcome measurements

  # Make Y NA if R==0
  Y[R==0] <- NA
  X[R==0] <- NA
  A[R==0] <- NA
  W <- X+rnorm(n)
  X <- as.data.frame(X)
  W <- as.data.frame(W)


  # too large of a cross-fitting value
  expect_error(check_entry_errors(Y,A,X,W,R,
                                  eem_ind,Rprobs,
                                  k=n+10))

  # X not in df format
  expect_error(check_entry_errors(Y,A,as.vector(X),
                                  W,R,
                                  eem_ind,Rprobs,
                                  k=1))

  # W not in df format
  expect_error(check_entry_errors(Y,A,X,
                                  as.vector(W),R,
                                  eem_ind,Rprobs,
                                  k=1))

  # Y not a numeric vector
  Ybad <- Y
  Ybad[1] <- 'a'
  expect_error(check_entry_errors(Ybad,A,X,W,R,
                                  eem_ind,Rprobs,
                                  k=1))

  # A not a numeric vector
  Abad <- A
  Abad[1] <- 'a'
  expect_error(check_entry_errors(Y,Abad,X,W,R,
                                  eem_ind,Rprobs,
                                  k=1))

  # catches naming conflicts (a few varnames arent allowed)
  for (bad_name in c('y','Y','A')) {
    colnames(X) <- bad_name

    expect_error(check_entry_errors(Y,A,X,W,R,
                                    eem_ind,Rprobs,
                                    k=1))
  }
  colnames(X) <- 'X'

  # User-supplied probs need to be valid probabilities
  expect_error(
    check_entry_errors(Y,A,X,W,R,
                                  eem_ind,Rprobs=rnorm(n),
                                  k=1)
  )

  # catches differences in length
  expect_error(check_entry_errors(rnorm(n-10),A,X,W,R,
                                  eem_ind,Rprobs,
                                  k=1))
  expect_error(check_entry_errors(Y,rnorm(n-10),X,W,R,
                                  eem_ind,Rprobs,
                                  k=1))

  # enforces eem_ind to be a logical
  expect_error(check_entry_errors(Y,A,X,W,R,
                                  'this should be a logical',Rprobs,
                                  k=1))



})

test_that('catches unspecified learners', {
  m_learners <- c('SL.glm','SL.mean')
  g_learners <- 'SL.glm'
  default_learners <- po_learners <- r_learners <- NULL

  expect_error(
    clean_learners(default_learners,
                   m_learners,g_learners,
                   r_learners,po_learners)
  )

})

test_that("create_folds returns correct structure for k=1", {
  folds <- create_folds(100, 1)
  expect_true(is.list(folds))
  expect_equal(folds$train, 1:100)
  expect_equal(folds$test, 1:100)
})

test_that("create_folds returns k folds for k>1", {
  folds <- create_folds(100, 5)
  expect_equal(length(folds), 5)

  # all indices should be covered exactly once across test sets
  all_test <- sort(unlist(lapply(folds, `[[`, "test")))
  expect_equal(all_test, 1:100)

  # train and test should be disjoint within each fold
  for (fold in folds) {
    expect_equal(length(intersect(fold$train, fold$test)), 0)
  }
})

test_that("check_binary identifies binary and non-binary", {
  expect_true(check_binary(c(0,1,0,1,1)))
  expect_true(check_binary(c(0, 0.5, 1))) # unit interval
  expect_false(check_binary(c(-1, 0, 1)))
  expect_false(check_binary(c(0, 1, 2)))
})

test_that("truncate_g clips extreme values", {
  x <- c(0.001, 0.5, 0.999)
  result <- suppressWarnings(truncate_g(x, cutoff = 0.025))
  expect_true(all(result >= 0.025))
  expect_true(all(result <= 0.975))
})

test_that("truncate_r clips extreme values", {
  x <- c(0.001, 0.5, 0.999)
  result <- suppressWarnings(truncate_r(x, cutoff = 0.01))
  expect_true(all(result >= 0.01))
  expect_true(all(result <= 0.99))
})

test_that("trim avoids exact 0 and 1", {
  x <- c(0, 0.5, 1)
  result <- trim(x)
  expect_true(result[1] > 0)
  expect_true(result[3] < 1)
  expect_equal(result[2], 0.5)
})

test_that("find_missing_pattern correctly identifies Z and U", {
  n <- 100
  X <- data.frame(X1 = rnorm(n))
  A <- rbinom(n,1,0.5)
  Y <- rnorm(n)
  Y[1:10] <- NA
  W <- X[,0]

  res <- find_missing_pattern(Y, A, X, W)

  expect_true("Y" %in% res$U | "y" %in% colnames(res$Z))
  expect_true("X1" %in% colnames(res$Z))
  expect_equal(sum(res$R), n - 10)
})

test_that("find_missing_pattern warns on very few complete cases", {
  n <- 1000
  X <- data.frame(X1 = rnorm(n))
  A <- rbinom(n,1,0.5)
  Y <- rnorm(n)
  Y[1:(n-5)] <- NA # only 5 complete cases
  W <- X[,0]

  expect_warning(find_missing_pattern(Y, A, X, W),
                 'Small number of complete cases')
})

test_that("clean_learners fills in defaults", {
  result <- clean_learners('SL.glm', NULL, NULL, NULL, NULL)
  expect_equal(result$m_learners, 'SL.glm')
  expect_equal(result$g_learners, 'SL.glm')
  expect_equal(result$r_learners, 'SL.glm')
  expect_equal(result$po_learners, 'SL.glm')
})

test_that("clean_learners preserves specific learners over defaults", {
  result <- clean_learners('SL.glm', 'SL.gam', NULL, NULL, NULL)
  expect_equal(result$m_learners, 'SL.gam')
  expect_equal(result$g_learners, 'SL.glm')
})

