context("DP")

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
alg=lol.project.dp
# Full-Rank Cases
n <- 100
d <- 6
K <- 2
data <- list(separable=lol.sims.rtrunk(n, d),
             unseparable=lol.sims.xor2(n, d))
# data piling isn't very good, so raise the threshold
p <- .5

test_that("DP works for full-rank", {
  run_data_test(data, alg=alg, p=p, sep=FALSE)
})

set.seed(123456)
# Low-Rank Cases
n <- 100
d <- 110
K <- 2
data <- list(separable=lol.sims.rtrunk(n, d),
             unseparable=lol.sims.xor2(n, d))

test_that("DP exhibits MDP for low-rank", {
  run_data_test(data, alg=alg, piled=TRUE)
})

n <- 100
d <- 3
test_that("DP fails with 3-class", {
  data <- lol.sims.mean_diff(n, d, K=d+1, md=3)
  expect_error(alg(data$X, data$Y, d-1))
})
