context("PCA")
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
alg=lol.project.pca
# Full-Rank Cases
n <- 100
d <- 6
K <- 2
data <- list(separable=lol.sims.rtrunk(n, d),
             unseparable=lol.sims.xor2(n, d))

test_that("PCA full-rank fails for r > d", {
  expect_error(alg(X=data$separable$X, Y=data$separable$Y, r=d+1))
})

test_that("PCA full-rank r == d", {
  r <- d
  run_data_test(data, alg=alg, r)
})

test_that("PCA full-rank r == K+1", {
  r <- K+1
  run_data_test(data, alg=alg, r)
})

test_that("PCA full-rank r == K", {
  r <- K
  run_data_test(data, alg=alg, r)
})

test_that("PCA full-rank r == 1", {
  r <- 1
  p=0.5
  run_data_test(data, alg=alg, r=r, p=p)
})


set.seed(12345)
# Low-Rank Cases
n <- 100
d <- 110
K <- 2
data <- list(separable=lol.sims.rtrunk(n, d),
             unseparable=lol.sims.xor2(n, d))
p=0.3

test_that("PCA low-rank fails for r > d", {
  expect_error(alg(X=data$separable$dat$X, Y=data$separable$dat$Y, r=d+1))
})

test_that("PCA low-rank r == K+1", {
  r <- K+1
  run_data_test(data, alg=alg, r)
})

test_that("PCA low-rank r == K", {
  r <- K
  run_data_test(data, alg=alg, r)
})

test_that("PCA low-rank r == 1", {
  r <- 1
  p=0.5
  run_data_test(data, alg=alg, r=r, p=p)
})

n <- 100
d <- 3
test_that("PCA works with multi-class", {
  data <- lol.sims.mean_diff(n, d, K=d+1, md=3)
  expect_error(alg(X=data$X, Y=data$Y, r=d-1), NA)
})

