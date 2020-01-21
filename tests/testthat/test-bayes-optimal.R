context("bayes-optimal")
library(MASS)

run_data_test <- function(data, alg, r=NULL, sep=TRUE, piled=FALSE, p=.05){
  result <- lapply(data, function(dat) {
    embed <- do.call(alg, list(X=dat$X, Y=dat$Y, r=r, mus=dat$mus,
                               Sigmas=dat$Sigmas, priors=dat$priors))
    if (piled) {
      classifier.alg = lol.classify.nearestCentroid
      classifier.return = NaN
    } else {
      classifier.alg = lda
      classifier.return = "class"
    }
    # run cross-validation with 10-fold validation
    class <- do.call(classifier.alg, list(embed$Xr[,1,drop=FALSE], dat$Y))
    pred <- predict(class, embed$Xr[,1,drop=FALSE])
    lhat <- sum(pred$class != dat$Y)/length(dat$Y)
    if (!isTRUE(piled)) {
      # check that the embedding works
      if (!is.null(r)) {
        expect_equal(dim(embed$Xr), c(n, r))
      }
    } else {
      # check that the embedding produces piling in the first dim
      for (ylab in unique(dat$Y)) {
        # round to ten-thousandths place for tolerance
        expect_equal(length(unique(round(embed$Xr[dat$Y == ylab, 1]*10000))), 1)
      }
    }
    return(lhat)
  })
  if (sep) {
    # use non-parametric test to evaluate performance on separable and unseparable example
    expect_lt(result$separable, result$unseparable)
  }
  return(result)
}

suppressWarnings(RNGversion("3.5.0"))
set.seed(123456)
alg=lol.project.bayes_optimal
# Full-Rank Cases
n <- 100
d <- 6
K <- 2
data <- list(separable=lol.sims.rtrunk(n, d),
             unseparable=lol.sims.xor2(n, d))
cutoff <- 0.35

test_that("Bayes Optimal full-rank", {
  run_data_test(data, alg=alg)
})


set.seed(12345)
# Low-Rank Cases
n <- 100
d <- 110
K <- 2
data <- list(separable=lol.sims.rtrunk(n, d),
             unseparable=lol.sims.xor2(n, d))

test_that("Bayes Optimal low-rank", {
  run_data_test(data, alg=alg)
})

n <- 100
d <- 3
test_that("Bayes Optimal works with multi-class", {
  data <- lol.sims.mean_diff(n, d, K=d+1, md=3)
  expect_error(alg(data$X, data$Y, mus=data$mus, Sigmas=data$Sigmas, priors=data$priors), NA)
})
