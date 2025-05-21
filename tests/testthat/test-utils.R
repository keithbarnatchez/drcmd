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
                                  k=n+10,nboot=10))

  # X not in df format
  expect_error(check_entry_errors(Y,A,as.vector(X),
                                  W,R,
                                  eem_ind,Rprobs,
                                  k=1,nboot=10))

  # W not in df format
  expect_error(check_entry_errors(Y,A,X,
                                  as.vector(W),R,
                                  eem_ind,Rprobs,
                                  k=1,nboot=10))

  # Y not a numeric vector
  Ybad <- Y
  Ybad[1] <- 'a'
  expect_error(check_entry_errors(Ybad,A,X,W,R,
                                  eem_ind,Rprobs,
                                  k=1,nboot=10))

  # A not a numeric vector
  Abad <- A
  Abad[1] <- 'a'
  expect_error(check_entry_errors(Y,Abad,X,W,R,
                                  eem_ind,Rprobs,
                                  k=1,nboot=10))

  # catches naming conflicts (a few varnames arent allowed)
  for (bad_name in c('y','Y','A')) {
    colnames(X) <- bad_name

    expect_error(check_entry_errors(Y,A,X,W,R,
                                    eem_ind,Rprobs,
                                    k=1,nboot=10))
  }
  colnames(X) <- 'X'

  # User-supplied probs need to be valid probabilities
  expect_error(
    check_entry_errors(Y,A,X,W,R,
                                  eem_ind,Rprobs=rnorm(n),
                                  k=1,nboot=10)
  )

  # catches differences in length
  expect_error(check_entry_errors(rnorm(n-10),A,X,W,R,
                                  eem_ind,Rprobs,
                                  k=1,nboot=10))
  expect_error(check_entry_errors(Y,rnorm(n-10),X,W,R,
                                  eem_ind,Rprobs,
                                  k=1,nboot=10))

  # enforces eem_ind to be a logical
  expect_error(check_entry_errors(Y,A,X,W,R,
                                  'this should be a logical',Rprobs,
                                  k=1,nboot=10))



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

