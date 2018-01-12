require(ggplot2)
require(lol)
require(reshape2)
require(Rmisc)
require(randomForest)
require(gridExtra)
require(latex2exp)

nsim <- 50

ds <- c(15, 50, 100, 500, 1000)
n <- 100
nd <- length(ds)

X <- list()
Y <- list()
for (i in 1:length(ds)) {
  d <- ds[i]
  Xe <- array(0, dim=c(nsim, n, d))
  Ye <- array(0, dim=c(nsim, n))
  for (j in 1:nsim) {
    res <- lol.sims.xor2(n, d, fall=100)
    Xe[j,,] <- res$X
    Ye[j,] <- res$Y
  }
  X[[i]] <- Xe
  Y[[i]] <- Ye
}
saveRDS(list(X=X, Y=Y), file='xor_sim.rds')

X <- list()
Y <- list()
for (i in 1:length(ds)) {
  d <- ds[i]
  Xe <- array(0, dim=c(nsim, n, d))
  Ye <- array(0, dim=c(nsim, n))
  for (j in 1:nsim) {
    res <- lol.sims.cigar(n, d)
    Xe[j,,] <- res$X
    Ye[j,] <- res$Y
  }
  X[[i]] <- Xe
  Y[[i]] <- Ye
}
saveRDS(list(X=X, Y=Y), file='cig_sim.rds')



X <- list()
Y <- list()
for (i in 1:length(ds)) {
  d <- ds[i]
  Xe <- array(0, dim=c(nsim, n, d))
  Ye <- array(0, dim=c(nsim, n))
  for (j in 1:nsim) {
    res <- lol.sims.rtrunk(n, d, r=TRUE)
    Xe[j,,] <- res$X
    Ye[j,] <- res$Y
  }
  X[[i]] <- Xe
  Y[[i]] <- Ye
}
saveRDS(list(X=X, Y=Y), file='rtr_sim.rds')



X <- list()
Y <- list()
for (i in 1:length(ds)) {
  d <- ds[i]
  Xe <- array(0, dim=c(nsim, n, d))
  Ye <- array(0, dim=c(nsim, n))
  for (j in 1:nsim) {
    res <- lol.sims.rtrunk(n, d, C=3, r=TRUE)
    Xe[j,,] <- res$X
    Ye[j,] <- res$Y
  }
  X[[i]] <- Xe
  Y[[i]] <- Ye
}
saveRDS(list(X=X, Y=Y), file='rtr3_sim.rds')


X <- list()
Y <- list()
for (i in 1:length(ds)) {
  d <- ds[i]
  Xe <- array(0, dim=c(nsim, n, d))
  Ye <- array(0, dim=c(nsim, n))
  for (j in 1:nsim) {
    res <- lol.sims.toep(n, d)
    Xe[j,,] <- res$X
    Ye[j,] <- res$Y
  }
  X[[i]] <- Xe
  Y[[i]] <- Ye
}
saveRDS(list(X=X, Y=Y), file='toe_sim.rds')


X <- list()
Y <- list()
for (i in 1:length(ds)) {
  d <- ds[i]
  Xe <- array(0, dim=c(nsim, n, d))
  Ye <- array(0, dim=c(nsim, n))
  for (j in 1:nsim) {
    res <- lol.sims.fat_tails(n, d, s0 = 1, r=TRUE)
    Xe[j,,] <- res$X
    Ye[j,] <- res$Y
  }
  X[[i]] <- Xe
  Y[[i]] <- Ye
}
saveRDS(list(X=X, Y=Y), file='ft_sim.rds')
