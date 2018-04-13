context("Nearest Centroid Classifier")

test_that("Nearest Centroid Finds correct Centers", {
  n <- 100
  d <- 6
  K <- 2
  nrep <- 10

  centroid.sim <- sapply(1:nrep, function(i) {
    dat=lol.sims.mean_diff(n, d, md=2)
    class <- lol.classify.nearestCentroid(dat$X, dat$Y)
    class$centroids
  })
  avg.centr.sim <- apply(centroid.sim, 1, FUN=mean)
  expect_equal(avg.centr.sim[1], 2, tolerance=0.1)
  expect_equal(avg.centr.sim[2:length(avg.centr.sim)], rep(0, length(avg.centr.sim)-1), tolerance=0.1)
})


test_that("Nearest Centroid Classifies Properly", {
  n <- 100
  d <- 2
  K <- 2
  nrep <- 10

  dat=lol.sims.mean_diff(n, d, md=3)
  class <- lol.classify.nearestCentroid(dat$X, dat$Y)
  Yhat <- predict(class, dat$X)
  lhat <- mean(Yhat != dat$Y)
  expect_lt(lhat, 0.1)
})
