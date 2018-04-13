context("Random Classifier")

set.seed(1234567)

test_that("Random Classifier Approximates Limit Error", {
  n <- 100
  d <- 6
  K <- 2
  nrep <- 10

  pr <- c(0.9, 0.1)
  lhat.sim <- sapply(1:nrep, function(i) {
    dat=lol.sims.mean_diff(n, d, priors=pr)
    class <- lol.classify.randomGuess(dat$X, dat$Y)
    expect_equal(class$priors, pr, tolerance=0.05)
    Yhat <- predict(class, dat$X)
    sum(Yhat == dat$Y)/length(Yhat)
  })

  expect_gte(mean(lhat.sim), 0.9^2)
})
