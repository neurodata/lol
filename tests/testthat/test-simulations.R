context("Simulations")

n <- 300
d <- 10

test_that("Fat Tails", {
  result <- lol.sims.fat_tails(n, d)
  expect_equal(dim(result$X)[1], n)
  expect_equal(dim(result$X)[2], d)
  expect_equal(length(result$Y), n)

  pr <- c(0.8, 0.2)
  result <- lol.sims.fat_tails(n, d, priors=pr)
  expect_equal(mean(result$Y == 1), pr[1], tolerance=0.05)
  expect_equal(mean(result$Y == 2), pr[2], tolerance=0.05)
  expect_equal(length(result$Y), n)

  expect_error(lol.sims.fat_tails(n, d, priors=c(0.9, 0.2)))
})



test_that("Mean Diff", {
  result <- lol.sims.mean_diff(n, d)
  expect_equal(dim(result$X)[1], n)
  expect_equal(dim(result$X)[2], d)
  expect_equal(length(result$Y), n)

  pr <- c(0.8, 0.2)
  result <- lol.sims.mean_diff(n, d, priors=pr)
  expect_equal(mean(result$Y == 1), pr[1], tolerance=0.05)
  expect_equal(mean(result$Y == 2), pr[2], tolerance=0.05)
  expect_equal(length(result$Y), n)

  expect_error(lol.sims.mean_diff(n, d, priors=c(0.9, 0.2)))
})


test_that("Toeplitz", {
  result <- lol.sims.toep(n, d)
  expect_equal(dim(result$X)[1], n)
  expect_equal(dim(result$X)[2], d)
  expect_equal(length(result$Y), n)

  pr <- c(0.8, 0.2)
  result <- lol.sims.toep(n, d, priors=pr)
  expect_equal(mean(result$Y == 1), pr[1], tolerance=0.05)
  expect_equal(mean(result$Y == 2), pr[2], tolerance=0.05)
  expect_equal(length(result$Y), n)

  expect_error(lol.sims.toep(n, d, priors=c(0.9, 0.2)))
})


test_that("Quadratic Discriminant Toeplitz", {
  result <- lol.sims.qdtoep(n, d)
  expect_equal(dim(result$X)[1], n)
  expect_equal(dim(result$X)[2], d)
  expect_equal(length(result$Y), n)

  pr <- c(0.8, 0.2)
  result <- lol.sims.qdtoep(n, d, priors=pr)
  expect_equal(mean(result$Y == 1), pr[1], tolerance=0.05)
  expect_equal(mean(result$Y == 2), pr[2], tolerance=0.05)
  expect_equal(length(result$Y), n)

  expect_error(lol.sims.qdtoep(n, d, priors=c(0.9, 0.2)))
})

test_that("Random Trunk", {
  result <- lol.sims.rtrunk(n, d)
  expect_equal(dim(result$X)[1], n)
  expect_equal(dim(result$X)[2], d)
  expect_equal(length(result$Y), n)

  pr <- c(0.8, 0.2)
  result <- lol.sims.rtrunk(n, d, priors=pr)
  expect_equal(mean(result$Y == 1), pr[1], tolerance=0.05)
  expect_equal(mean(result$Y == 2), pr[2], tolerance=0.05)
  expect_equal(length(result$Y), n)

  expect_error(lol.sims.rtrunk(n, d, priors=c(0.9, 0.2)))

  result <- lol.sims.rtrunk(n, d, K=3)
  expect_equal(length(unique(result$Y)), 3)
  expect_equal(dim(result$X)[1], n)
  expect_equal(dim(result$X)[2], d)
  expect_equal(length(result$Y), n)
})

test_that("Cigar", {
  result <- lol.sims.cigar(n, d)
  expect_equal(dim(result$X)[1], n)
  expect_equal(dim(result$X)[2], d)
  expect_equal(length(result$Y), n)

  pr <- c(0.8, 0.2)
  result <- lol.sims.cigar(n, d, priors=pr)
  expect_equal(mean(result$Y == 1), pr[1], tolerance=0.05)
  expect_equal(mean(result$Y == 2), pr[2], tolerance=0.05)
  expect_equal(length(result$Y), n)

  expect_error(lol.sims.cigar(n, d, priors=c(0.9, 0.2)))
})

test_that("XOR", {
  result <- lol.sims.xor2(n, d)
  expect_equal(dim(result$X)[1], n)
  expect_equal(dim(result$X)[2], d)
  expect_equal(length(result$Y), n)

  pr <- c(0.8, 0.2)
  result <- lol.sims.xor2(n, d, priors=pr)
  expect_equal(mean(result$Y == 1), pr[1], tolerance=0.05)
  expect_equal(mean(result$Y == 2), pr[2], tolerance=0.05)
  expect_equal(length(result$Y), n)

  expect_error(lol.sims.xor2(n, d, priors=c(0.9, 0.2)))
})

test_that("Rotation", {
  pr <- c(0.8, 0.2)
  result <- lol.sims.mean_diff(n, d, md=3, priors=pr, rotate=TRUE)
  expect_equal(mean(result$Y == 1), pr[1], tolerance=0.05)
  expect_equal(mean(result$Y == 2), pr[2], tolerance=0.05)
  expect_equal(length(result$Y), n)
  expect_equal(dim(result$X),  c(n, d))

  result <- lol.sims.cigar(n, d, rotate=TRUE)
  expect_equal(dim(result$X)[1], n)
  expect_equal(dim(result$X)[2], d)
  expect_equal(length(result$Y), n)

  result <- lol.sims.fat_tails(n, d, rotate=TRUE)
  expect_equal(dim(result$X)[1], n)
  expect_equal(dim(result$X)[2], d)
  expect_equal(length(result$Y), n)

  result <- lol.sims.toep(n, d, rotate=TRUE)
  expect_equal(dim(result$X)[1], n)
  expect_equal(dim(result$X)[2], d)
  expect_equal(length(result$Y), n)

  result <- lol.sims.qdtoep(n, d, rotate=TRUE)
  expect_equal(dim(result$X)[1], n)
  expect_equal(dim(result$X)[2], d)
  expect_equal(length(result$Y), n)

  result <- lol.sims.rtrunk(n, d, rotate=TRUE)
  expect_equal(dim(result$X)[1], n)
  expect_equal(dim(result$X)[2], d)
  expect_equal(length(result$Y), n)
})
