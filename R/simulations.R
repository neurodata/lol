#' Fat Tails Simulation
#'
#' A function for simulating from 2 classes, one of which has fatter tails.
#'
#' @import abind
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param r=FALSE whether to apply a random rotation.
#' @param f=15 the fatness scaling of the tail.
#' @param s0=10 the number of dimensions with a difference in the means. s0 should be < d.
#' @param rho=0.2 the scaling of the covariance terms, should be < 1.
#' @return X the data as a [n, d] matrix.
#' @return Y the labels as a [d] array.
#' @author Eric Bridgeford, adapted from Joshua Vogelstein
#' @examples
#' library(fselect)
#' data <- fs.sims.fat_tails(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' @export
fs.sims.fat_tails <- function(n, d, r=FALSE, f=15, s0=10, rho=0.2) {
  if  (D1 > d) {
    stop(sprintf("s0 = %d, d=%d. s0 should be > d.", s0, d))
  }
  mu0 <- array(0, dim=c(d, 1))
  mu1 <- c(array(0, dim=c(s0)), array(1, dim=c(d - s0)))
  mus <- abind::abind(mu0, mu1, along=2)

  S <- rho*array(1, dim=c(d, d))
  diag(S) <- 1

  S <- abind::abind(S, 15*S, along=3)

  if (r) {
    res <- fs.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
  }

  # simulate from GMM
  sim <- fs.sims.sim_gmm(mus, S, n)
  return(list(X=sim$X, Y=sim$Y))
}

#' Toeplitz Simulation
#'
#' A function for simulating data in which the covariance is a non-symmetric toeplitz matrix.
#'
#' @import abind
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param D1=10 the dimensionality for the non-equal covariance terms.
#' @param r=FALSE whether to apply a random rotation.
#' @param b=0.4 a scaling parameter for the means.
#' @param rho=0.5 the scaling of the covariance terms, should be < 1.
#' @return X [n, d] the data as a matrix.
#' @return Y [n] the labels as a array.
#' @author Eric Bridgeford, adapted from Joshua Vogelstein
#' @examples
#' library(fselect)
#' data <- fs.sims.toep(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' @export
fs.sims.toep <- function(n, d, D1=10, r=FALSE, b=0.4, rho=0.5) {
  c <- rho^(0:(D1 - 1))
  A <- toeplitz(c)
  K1 <- sum(A)

  c <- rho^(0:(d-1))
  A <- toeplitz(c)
  K <- sum(A)

  mudelt <- (K1 * b^2/K)^0.5/2

  mu0 <- array(mudelt*c(1, -1), dim=c(d))
  mus <- abind::abind(mu0, -mu0, along=2)
  S <- abind::abind(A, A, along=3)

  if (r) {
    res <- fs.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
  }
  # simulate from GMM
  sim <- fs.sims.sim_gmm(mus, S, n)
  return(list(X=sim$X, Y=sim$Y))
}

#' Quadratic Discriminant Toeplitz Simulation
#'
#' A function for simulating data generalizing the Toeplitz setting, where each class has a different covariance matrix. This results in a Quadratic Discriminant.
#'
#' @import abind
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param D1=10 the dimensionality for the non-equal covariance terms.
#' @param r=FALSE whether to apply a random rotation.
#' @param b=0.4 a scaling parameter for the means.
#' @param rho=0.5 the scaling of the covariance terms, should be < 1.
#' @return X [n, d] the data as a matrix.
#' @return Y [n] the labels as a array.
#' @author Eric Bridgeford, adapted from Joshua Vogelstein
#' @examples
#' library(fselect)
#' data <- fs.sims.qdtoep(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' @export
fs.sims.qdtoep <- function(n, d, D1=10, r=FALSE, b=0.4, rho=0.5) {
  tR <- rho^(0:(D1 - 1))
  A <- toeplitz(tR)
  K1 <- sum(A)

  tR <- rho^(0:(d-1))
  A0 <- toeplitz(tR)
  K <- sum(A0)

  mudelt <- (K1 * b^2/K)^0.5/2

  mu0 <- array(mudelt*c(1, -1), dim=c(d))
  Q <- fs.sims.rotation(d)
  mu1 <- -Q %*% (mu0 + 0.1)
  mus <- abind::abind(mu0, mu1, along=2)
  A1 <- Q %*% A0 %*% t(Q)
  S <- abind::abind(A0, A1, along=3)

  if (r) {
    res <- fs.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
  }
  # simulate from GMM
  sim <- fs.sims.sim_gmm(mus, S, n)
  return(list(X=sim$X, Y=sim$Y))
}

#' Random Trunk
#'
#' A simulation for the random trunk experiment.
#' @import abind
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param b=4 scalar for mu scaling.
#' @param r=FALSE whether to apply a random rotation.
#' @param C=2 number of classes, should be <4.
#' @return X [n, d] the data as a matrix.
#' @return Y [n] the labels as a array.
#' @author Eric Bridgeford, adapted from Joshua Vogelstein
#' @examples
#' library(fselect)
#' data <- fs.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' @export
fs.sims.rtrunk <- function(n, d, b=4, r=FALSE, C=2) {
  mu1 <- b/sqrt(0:(d-1)*2 + 1)
  if (C == 2) {
    mus <- abind::abind(mu1, -mu1, along=2)
  } else if (C == 3) {
    mus <- abind::abind(mu1, 0*mu1, -mu1, along=2)
  } else if (C == 4) {
    mus <- abind::abind(mu1, b/(seq(from=d, to=1, by=-1)), b/(seq(from=1, to=d, by=1)), -mu1, along=2)
  }
  S <- diag(d)
  diag(S) <- 100/sqrt(seq(from=d, to=1, by=-1))

  S <- array(unlist(replicate(C, S, simplify=FALSE)), dim=c(d, d, C))

  if (r) {
    res <- fs.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
  }
  # simulate from GMM
  sim <- fs.sims.sim_gmm(mus, S, n)
  return(list(X=sim$X, Y=sim$Y))
}

#' Stacked Cigar
#'
#' A simulation for the stacked cigar experiment.
#' @import abind
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param a=0.15 scalar for all of the mu1 but 2nd dimension.
#' @param b=4 scalar for 2nd dimension value of mu2.
#' @param r=FALSE whether to apply a random rotation.
#' @return X [n, d] the data as a matrix.
#' @return Y [n] the labels as a array.
#' @author Eric Bridgeford, adapted from Joshua Vogelstein
#' @examples
#' library(fselect)
#' data <- fs.sims.cigar(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' @export
fs.sims.cigar <- function(n, d, a=0.15, b=4, r=FALSE) {
  mu1 <- array(a, dim=c(d))
  mu1[2] <- b
  mus <- cbind(array(0, dim=c(d)), mu1)

  S <- diag(d)
  S[2,2] <- b

  S <- abind::abind(diag(d), S, along=3)

  if (r) {
    res <- fs.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
  }
  # simulate from GMM
  sim <- fs.sims.sim_gmm(mus, S, n)
  return(list(X=sim$X, Y=sim$Y))
}


#' Xor Problem
#'
#' A function to simulate from the 2-class xor problem.
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param fall=100 the sigma for the covariance structuring.
#' @return X [n, d] the data as a matrix.
#' @return Y [n] the labels as a array.
#' @author Eric Bridgeford, adapted from Joshua Vogelstein
#' @examples
#' library(fselect)
#' data <- fs.sims.xor2(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' @export
fs.sims.xor2 <- function(n, d, fall=100) {
  n1 <- ceiling(n/2)
  n2 <- floor(n/2)
  # first simulation set
  mus <- abind::abind(array(0, dim=c(d)), array(c(1, 0), dim=c(d)), along=2)
  S <- sqrt(d/fall)*diag(d)
  S <- abind::abind(S, S, along=3)

  # simulate from GMM for first set of training examples
  sim1 <- fs.sims.sim_gmm(mus, S, n1)

  # second simulation set
  mus <- abind::abind(array(1, dim=c(d)), array(c(0, 1), dim=c(d)), along=2)

  # simulate from GMM for second set of training examples
  sim2 <- fs.sims.sim_gmm(mus, S, n2)

  X <- abind::abind(sim1$X, sim2$X, along=1)
  Y <- abind::abind(sim1$Y, sim2$Y, along=1)

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
  ylabs <- as.vector(sort(unique(labs)))
  res <- sapply(ylabs, function(y) mvrnorm(n=sum(labs == y), mus[,y], Sigmas[,,y]), USE.NAMES=TRUE, simplify=FALSE)
  X <- array(0, dim=c(n, dim(Sigmas)[1]))
  for (y in ylabs) {
    X[labs == y,] <- res[[y]]
  }
  return(list(X=X, Y=labs))
}

#' Sample Random Rotation
#'
#' A helper function for estimating a random rotation matrix.
#' @author Eric Bridgeford
fs.sims.rotation <- function(d) {
  Q <- Matrix::qr.Q(Matrix::qr(array(rnorm(d*d), dim=c(d, d))))
  if (det(Q) < -.99) {
    Q[,1] <- -Q[,1]
  }
  return(Q)
}

#' Random Rotation
#'
#' A helper function for applying a random rotation to gaussian parameter set.
#' @author Eric Bridgeford
fs.sims.random_rotate <- function(mus, Sigmas) {
  dimm <- dim(mus)
  C <- dimm[2]
  d <- dim(mus)[1]
  Q <- fs.sims.rotation(d)

  for (i in 1:C) {
    mus[,i] <- Q %*% mus[,i,drop=FALSE]
    Sigmas[,,i] <- Q %*% Sigmas[,,i] %*% t(Q)
  }
  return(list(mus=mus, S=Sigmas))
}
