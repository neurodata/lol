context("Cross Validation")

set.seed(123456)
n <- 100
d <- 20
xv=5
# labels and Xs don't have proper dimensions
dat <- lol.sims.mean_diff(n, d, md=2)
sets <- lol.xval.split(dat$X, dat$Y, k=xv)
alg <- lol.project.lol

test_that('Failure Casees', {
  r <- 2
  expect_error(lol.xval.eval(dat$X, c(dat$Y, dat$Y), r, lol.project.lol, sets=sets))
  expect_error(lol.xval.optimal_dimselect(dat$X, c(dat$Y, dat$Y), c(2, 20), lol.project.lol, sets=sets))

})

test_that('Cross Validation Returns Lhat', {
  r <- 2
  result <- lol.xval.eval(dat$X, dat$Y, r, alg=lol.project.lol, k=xv)
  expect_equal(length(result$lhats), xv)

  # works with any classifier
  result <- lol.xval.eval(dat$X, dat$Y, r, alg=lol.project.lol, k=xv,
                          classifier=lol.classify.nearestCentroid, classifier.return=NaN)
  expect_equal(length(result$lhats), xv)
})

test_that('Optimal Embedding Dimensions', {
  r <- 2
  result.r <- lol.xval.eval(dat$X, dat$Y, r, alg=lol.project.lol, sets=sets, k=xv)
  rs <- c(2, 20)
  result.rs <- lol.xval.optimal_dimselect(dat$X, dat$Y, rs, lol.project.lol, sets=sets, k=xv)
  expect_equal(r, result.rs$optimal.r)
  expect_equal(as.numeric(result.rs$folds.data$lhat[result.rs$folds.data$r == r]),
               as.numeric(result.r$lhats))

  # works when not structured
  result.rs.notstruct <- lol.xval.optimal_dimselect(dat$X, dat$Y, rs, lol.project.lol, sets=sets, k=xv, alg.structured = FALSE)
  expect_equal(as.numeric(result.rs$folds.data$lhat), as.numeric(result.rs.notstruct$folds.data$lhat))

  # no specified set
  result.rs.noset <- lol.xval.optimal_dimselect(dat$X, dat$Y, rs, lol.project.lol, k=xv, alg.structured = FALSE)
  expect_equal(result.rs$optimal.lhat, result.rs.noset$optimal.lhat, tolerance=0.05)

  # works with any classifier
  result <- lol.xval.optimal_dimselect(dat$X, dat$Y, rs, alg=lol.project.lol, k=xv,
                                       classifier=lol.classify.nearestCentroid, classifier.return=NaN)
  expect_equal(length(result$folds.data$lhat), xv*length(rs))
})

test_that('Kfold Setup', {
  sets <- lol.xval.split(dat$X, dat$Y, k='loo')
  expect_true(length(sets) == n)

  sets <- lol.xval.split(dat$X, dat$Y, k=20)
  expect_true(length(sets) == 20)

  expect_error(lol.xval.split(dat$X, dat$Y, k=NULL))
})
