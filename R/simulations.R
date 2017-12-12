#' Fat Tails Simulation
#'
#' A function for simulating from 2 classes, one of which has fatter tails.
#'
#' @import abind
#' @param D the dimensionality of the simulated data.
#' @param n the number of samples of the simulated data.
#' @param r=FALSE whether to apply a random rotation.
#' @param f=15 the fatness scaling of the tail.
#' @param s0=10 the number of dimensions with a difference in the means.
#' @param rho=0.2 the scaling of the covariance terms, should be < 1.
#' @return X the data as a [n, d] matrix.
#' @return Y the labels as a [d] array.
#' @author Eric Bridgeford, adapted from Joshua Vogelstein
#' @export
fs.sims.fat_tails <- function(D, n, r=FALSE, f=15, s0=10, t=0.8, rho=0.2) {
  mu0 <- array(0, dim=c(D, 1))
  mu1 <- c(array(1, dim=c(s0)), array(0, dim=c(D - s0)))
  mus <- abind(mu0, mu1, along=2)

  S <- rho*array(1, dim=c(D, D))
  diag(S) <- 1
  if (r) {
    res <- fs.sims.random_rotate(mu, S)
    mus <- res$mu
    S <- res$S
  }

  S <- abind(S, 15*S, along=3)

  # simulate from GMM
  sim <- fs.sims.sim_gmm(mus, S, n)
  return(list(X=sim$X, Y=sim$Y))
}

#' Toeplitz Simulation
#'
#' A function for simulating data in which the covariance is a non-symmetric toeplitz matrix.
#'
#' @import abind
#' @param D the dimensionality of the simulated data.
#' @param n the number of samples of the simulated data.
#' @param D1=10 the dimensionality for the non-equal covariance terms
#' @param r=FALSE whether to apply a random rotation.
#' @param rho=0.5 the scaling of the covariance terms, should be < 1.
#' @return X [n, d] the data as a matrix.
#' @return Y [n] the labels as a array.
#' @author Eric Bridgeford, adapted from Joshua Vogelstein
#' @export
fs.sims.toep <- function(D, n, D1=10, r=FALSE, b=0.4, rho=0.5) {
  c <- rho^(0:(D1 - 1))
  A <- toeplitz(c)
  K1 <- sum(A)

  c <- rho^(0:(D-1))
  A <- toeplitz(c)
  K <- sum(A)

  mudelt <- (K1 * b^2/K)^0.5/2

  mu0 <- array(mudelt*c(1, -1), dim=c(D))
  mus <- abind(mu0, -mu0, along=2)
  S <- abind(A, A, along=3)

  # simulate from GMM
  sim <- fs.sims.sim_gmm(mus, S, n)
  return(list(X=sim$X, Y=sim$Y))
}

#' Random Trunk
#'
#' A simulation for the random trunk experiment.
#' @import abind
#' @param D the dimensionality of the simulated data.
#' @param n the number of samples of the simulated data.
#' @param b=4 scalar for mu scaling.
#' @param r=FALSE whether to apply a random rotation.
#' @param C=2 number of classes, should be <4.
#' @return X [n, d] the data as a matrix.
#' @return Y [n] the labels as a array.
#' @author Eric Bridgeford, adapted from Joshua Vogelstein
#' @export
fs.sims.rtrunk <- function(D, n, b=4, r=FALSE, C=2) {
  mu1 <- b/sqrt(0:(D-1)*2 + 1)
  if (C == 2) {
    mus <- abind(mu1, -mu1, along=2)
  } else if (C == 3) {
    mus <- abind(mu1, 0*mu1, -mu1, along=2)
  } else if (C == 4) {
    mus <- abind(mu1, b/(seq(from=D, to=1, by=-1)), b/(seq(from=1, to=D, by=1)), -mu1, along=3)
  }
  S <- diag(D)
  diag(S) <- 100/sqrt(seq(from=D, to=1, by=-1))

  if (r) {
    res <- fs.sims.random_rotate(mus, S)
    mus <- res$mu
    S <- res$S
  }
  S <- array(unlist(replicate(C, S, simplify=FALSE)), dim=c(D, D, C))

  # simulate from GMM
  sim <- fs.sims.sim_gmm(mus, S, n)
  return(list(X=sim$X, Y=sim$Y))
}

#' Xor Problem
#'
#' A function to simulate from the 2-class xor problem.
#' @param D the dimensionality of the simulated data.
#' @param n the number of samples of the simulated data.
#' @return X [n, d] the data as a matrix.
#' @return Y [n] the labels as a array.
#' @author Eric Bridgeford, adapted from Joshua Vogelstein
#' @export
fs.sims.xor2 <- function(D, n) {
  n1 <- ceiling(n/2)
  n2 <- floor(n/2)
  # first simulation set
  mus <- abind(array(0, dim=c(D)), array(c(1, 0), dim=c(D)), along=2)
  S <- sqrt(D/4)*diag(D)
  S <- abind(S, S, along=3)

  # simulate from GMM for first set of training examples
  sim1 <- fs.sims.sim_gmm(mus, S, n1)

  # second simulation set
  mus <- abind(array(1, dim=c(D)), array(c(0, 1), dim=c(D)), along=2)
  S <- sqrt(D/4)*diag(D)
  S <- abind(S, S, along=3)

  # simulate from GMM for second set of training examples
  sim2 <- fs.sims.sim_gmm(mus, S, n2)

  X <- abind(sim1$X, sim2$X, along=1)
  Y <- abind(sim1$Y, sim2$Y, along=1)

  reorder <- sample(n)
  return(list(X=X[reorder,], Y=Y[reorder]))
}

#' GMM Simulate
#'
#' A helper function for simulating from Gaussian Mixture.
#' @param mus [d, C] the mus for each class.
#' @param Sigmas [d,d,C] the Sigmas for each class.
#' @param n the of examples.
#' @return X [n, d] the simulated data.
#' @return Y [n] the labels for each data point.
#' @author Eric Bridgeford
#' @import MASS
fs.sims.sim_gmm <- function(mus, Sigmas, n) {
  C <- dim(mus)[2]
  labs <- sample(1:C, size=n, replace=TRUE)
  X <- sapply(labs, function(lab) mvrnorm(n=1, mus[,lab], Sigmas[,,lab]), simplify=FALSE)
  return(list(X=matrix(unlist(X), nrow=length(X), byrow=FALSE), Y=labs))
}

#' Random Rotation
#'
#' A helper function for applying a random rotation to gaussian parameter set.
#' @author Eric Bridgeford, adapted from Joshua Vogelstein
fs.sims.random_rotate <- function(mu, S) {
  stop("Not Implemented!")
}
