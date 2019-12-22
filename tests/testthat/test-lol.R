context("LOL")

library(MASS)

run_data_test <- function(data, alg, r=NULL, sep=TRUE, piled=FALSE, p=.05){
  result <- lapply(data, function(dat) {
    if (piled) {
      classifier.alg = lol.classify.nearestCentroid
      classifier.return = NaN
    } else {
      classifier.alg = lda
      classifier.return = "class"
    }
    # run cross-validation with 10-fold validation
    result <- expect_error(lol.xval.eval(dat$X, dat$Y, r, alg=alg, classifier=classifier.alg,
                                         classifier.return=classifier.return, k=10), NA)
    if (!isTRUE(piled)) {
      # check that the embedding works
      if (!is.null(r)) {
        expect_equal(dim(result$model$Xr), c(n, r))
      }
    } else {
      # check that the embedding produces piling in the first dim
      for (ylab in unique(dat$Y)) {
        # round to ten-thousandths place for tolerance
        expect_equal(length(unique(round(result$model$Xr[dat$Y == ylab, 1]*10000))), 1)
      }
    }
    return(result$lhats)
  })
  if (sep) {
    # use non-parametric test to evaluate performance on separable and unseparable example
    expect_lt(wilcox.test(result$separable, result$unseparable, alternative="less", exact=FALSE)$p.value, p)
  }
  return(result)
}

suppressWarnings(RNGversion("3.5.0"))
set.seed(123456)
alg=lol.project.lol
# Full-Rank Cases
n <- 100
d <- 6
K <- 2
data <- list(separable=lol.sims.rtrunk(n, d),
             unseparable=lol.sims.xor2(n, d))

test_that("LOL full-rank fails for r > d", {
  expect_error(alg(data$separable$X, data$separable$Y, d+1))
})

test_that("LOL full-rank r == d", {
  r <- d
  run_data_test(data, alg=alg, r)
})

test_that("LOL full-rank r == K+1", {
  r <- K+1
  run_data_test(data, alg=alg, r)
})

test_that("LOL full-rank r == K", {
  r <- K
  run_data_test(data, alg=alg, r)
})

test_that("LOL full-rank r == 1", {
  r <- 1
  run_data_test(data, alg=alg, r=r)
})


set.seed(12345)
# Low-Rank Cases
n <- 100
d <- 110
K <- 2
data <- list(separable=lol.sims.rtrunk(n, d),
             unseparable=lol.sims.xor2(n, d))

test_that("LOL low-rank fails for r > d", {
  expect_error(alg(data$separable$ data$separable$Y, d+1))
})

test_that("LOL low-rank r == K+1", {
  r <- K+1
  run_data_test(data, alg=alg, r)
})

test_that("LOL low-rank r == K", {
  r <- K
  run_data_test(data, alg=alg, r)
})

test_that("LOL low-rank r == 1", {
  r <- 1
  run_data_test(data, alg=alg, r=r)
})

n <- 100
d <- 3
test_that("LOL works with multi-class", {
  data <- lol.sims.mean_diff(n, d, K=d+1, md=3)
  expect_error(alg(data$X, data$Y, d-1), NA)
})
