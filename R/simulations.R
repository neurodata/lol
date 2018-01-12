#' Fat Tails Simulation
#'
#' A function for simulating from 2 classes, one of which has fatter tails.
#'
#' @import abind
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param rotate=FALSE whether to apply a random rotation to the mean and covariance. With random rotataion matrix Q, mu = Q*mu, and S = Q*S*Q.
#' @param priors=NULL the priors for each class. If NULL, class priors are all equal. If not null, should be |priors| = K, a length K vector for K classes.
#' @param f=15 the fatness scaling of the tail. S2 = f*S1, where S1_{ij} = rho if i != j, and 1 if i == j.
#' @param s0=10 the number of dimensions with a difference in the means. s0 should be < d.
#' @param rho=0.2 the scaling of the off-diagonal covariance terms, should be < 1.
#' @return X the data as a [n, d] matrix.
#' @return Y the labels as a [d] array.
#' @return mus [d, K] the per-class means.
#' @return Sigmas [d, d, K] the per-class covariance matrices.
#' @author Eric Bridgeford
#' @examples
#' library(lol)
#' data <- lol.sims.fat_tails(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' @export
lol.sims.fat_tails <- function(n, d, rotate=FALSE, f=15, s0=10, rho=0.2, priors=NULL) {
  if  (s0 > d) {
    stop(sprintf("s0 = %d, d=%d. s0 should be < d.", s0, d))
  }
  mu0 <- array(0, dim=c(d, 1))
  mu1 <- c(array(0, dim=c(s0)), array(1, dim=c(d - s0)))
  mus <- abind::abind(mu0, mu1, along=2)

  S <- array(rho, dim=c(d, d))
  diag(S) <- 1

  S <- abind::abind(S, 15*S, along=3)

  if (rotate) {
    res <- lol.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
  }

  # simulate from GMM
  sim <- lol.sims.sim_gmm(mus, S, n, priors)
  return(list(X=sim$X, Y=sim$Y, mus=mus, Sigmas=S))
}

#' Mean Difference Simulation
#'
#' A function for simulating data in which a difference in the means is present only in a subset of dimensions, and equal covariance.
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param rotate=FALSE whether to apply a random rotation to the mean and covariance. With random rotataion matrix Q, mu = Q*mu, and S = Q*S*Q.
#' @param priors=NULL the priors for each class. If NULL, class priors are all equal. If not null, should be |priors| = K, a length K vector for K classes.
#' @param K=2 the number of classes.
#' @param md=1 the magnitude of the difference in the means in the specified subset of dimensions.
#' @param subset=c(1) the dimensions to have a difference in the means. Defaults to only the first dimension. max(subset) < d.
#' @param offdiag=0 the off-diagonal elements of the covariance matrix. Should be < 1. S_{ij} = offdiag if i != j, or 1 if i == j.
#' @param s=1 the scaling parameter of the covariance matrix. S_{ij} = scaling*1 if i == j, or scaling*offdiag if i != j.
#' @return X [n, d] the data as a matrix.
#' @return Y [n] the labels as a array.
#' @return mus [d, K] the per-class means.
#' @return Sigmas [d, d, K] the per-class covariance matrices.
lol.sims.mean_diff <- function(n, d, rotate=FALSE, priors=NULL, K=2, md=1, subset=c(1), offdiag=0, s=1) {
  if (max(subset) > d) {
    stop(sprintf("Specified a difference in dimension %d; maximum should be %d.", max(subset), d))
  }
  mus <- array(0, dim=c(d, K))
  for (i in 1:K) {
    mus[subset] <- i*md  # scale the difference in the means accordingly
  }
  S <- array(offdiag, dim=c(d, d))
  diag(S) <- 1  # identity variance
  S <- s*S  # scale accordingly

  S <- array(unlist(replicate(K, S, simplify=FALSE)), dim=c(d, d, K))

  if (rotate) {
    res <- lol.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
  }
  # simulate from GMM
  sim <- lol.sims.sim_gmm(mus, S, n, priors)
  return(list(X=sim$X, Y=sim$Y, mus=mus, Sigmas=S))
}

#' Toeplitz Simulation
#'
#' A function for simulating data in which the covariance is a non-symmetric toeplitz matrix.
#'
#' @import abind
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param D1=10 the dimensionality for the non-equal covariance terms.
#' @param rotate=FALSE whether to apply a random rotation to the mean and covariance. With random rotataion matrix Q, mu = Q*mu, and S = Q*S*Q.
#' @param priors=NULL the priors for each class. If NULL, class priors are all equal. If not null, should be |priors| = K, a length K vector for K classes.
#' @param b=0.4 a scaling parameter for the means.
#' @param rho=0.5 the scaling of the covariance terms, should be < 1.
#' @return X [n, d] the data as a matrix.
#' @return Y [n] the labels as a array.
#' @return mus [d, K] the per-class means.
#' @return Sigmas [d, d, K] the per-class covariance matrices.
#' @author Eric Bridgeford
#' @examples
#' library(lol)
#' data <- lol.sims.toep(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' @export
lol.sims.toep <- function(n, d, rotate=FALSE, priors=NULL, D1=10, b=0.4, rho=0.5) {
  if (rho >= 1) {
    stop(sprintf("rho should be < 1; user specified %.3f", rho))
  }
  cT <- rho^(0:(D1 - 1))
  A <- toeplitz(cT)
  K1 <- sum(A)

  cT <- rho^(0:(d-1))
  A <- toeplitz(cT)
  K <- sum(A)

  mudelt <- (K1 * b^2/K)^0.5/2

  mu0 <- array(mudelt*c(1, -1), dim=c(d))
  mus <- abind::abind(mu0, -mu0, along=2)
  S <- abind::abind(A, A, along=3)

  if (rotate) {
    res <- lol.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
  }
  # simulate from GMM
  sim <- lol.sims.sim_gmm(mus, S, n, priors)
  return(list(X=sim$X, Y=sim$Y, mus=mus, Sigmas=S))
}

#' Quadratic Discriminant Toeplitz Simulation
#'
#' A function for simulating data generalizing the Toeplitz setting, where each class has a different covariance matrix. This results in a Quadratic Discriminant.
#'
#' @import abind
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param D1=10 the dimensionality for the non-equal covariance terms.
#' @param rotate=FALSE whether to apply a random rotation to the mean and covariance. With random rotataion matrix Q, mu = Q*mu, and S = Q*S*Q.
#' @param priors=NULL the priors for each class. If NULL, class priors are all equal. If not null, should be |priors| = K, a length K vector for K classes.
#' @param b=0.4 a scaling parameter for the means.
#' @param rho=0.5 the scaling of the covariance terms, should be < 1.
#' @return X [n, d] the data as a matrix.
#' @return Y [n] the labels as a array.
#' @return mus [d, K] the per-class means.
#' @return Sigmas [d, d, K] the per-class covariance matrices.
#' @author Eric Bridgeford
#' @examples
#' library(lol)
#' data <- lol.sims.qdtoep(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' @export
lol.sims.qdtoep <- function(n, d, rotate=FALSE, priors=NULL, D1=10, b=0.4, rho=0.5) {
  if (rho >= 1) {
    stop(sprintf("rho should be < 1; user specified %.3f", rho))
  }
  if (D1 > d) {
    stop(sprintf("User has specified more dimensions for non-equal cov terms. D1 is %d, yet d is %d", D1, d.))
  }
  tR <- rho^(0:(D1 - 1))
  A <- toeplitz(tR)
  K1 <- sum(A)

  tR <- rho^(0:(d-1))
  A0 <- toeplitz(tR)
  K <- sum(A0)

  mudelt <- (K1 * b^2/K)^0.5/2

  mu0 <- array(mudelt*c(1, -1), dim=c(d))
  Q <- lol.sims.rotation(d)
  mu1 <- -Q %*% (mu0 + 0.1)
  mus <- abind::abind(mu0, mu1, along=2)
  A1 <- Q %*% A0 %*% t(Q)
  S <- abind::abind(A0, A1, along=3)

  if (rotate) {
    res <- lol.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
  }
  # simulate from GMM
  sim <- lol.sims.sim_gmm(mus, S, n, priors)
  return(list(X=sim$X, Y=sim$Y, mus=mus, Sigmas=S))
}

#' Random Trunk
#'
#' A simulation for the random trunk experiment.
#' @import abind
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param rotate=FALSE whether to apply a random rotation to the mean and covariance. With random rotataion matrix Q, mu = Q*mu, and S = Q*S*Q.
#' @param priors=NULL the priors for each class. If NULL, class priors are all equal. If not null, should be |priors| = K, a length K vector for K classes.
#' @param b=4 scalar for mu scaling.
#' @param K=2 number of classes, should be <4.
#' @return X [n, d] the data as a matrix.
#' @return Y [n] the labels as a array.
#' @return mus [d, K] the per-class means.
#' @return Sigmas [d, d, K] the per-class covariance matrices.
#' @author Eric Bridgeford
#' @examples
#' library(lol)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' @export
lol.sims.rtrunk <- function(n, d, rotate=FALSE, priors=NULL, b=4, K=2) {
  mu1 <- b/sqrt(0:(d-1)*2 + 1)
  if (K == 2) {
    mus <- abind::abind(mu1, -mu1, along=2)
  } else if (K == 3) {
    mus <- abind::abind(mu1, 0*mu1, -mu1, along=2)
  } else if (K == 4) {
    mus <- abind::abind(mu1, b/(seq(from=d, to=1, by=-1)), b/(seq(from=1, to=d, by=1)), -mu1, along=2)
  }
  S <- diag(d)
  diag(S) <- 100/sqrt(seq(from=d, to=1, by=-1))

  S <- array(unlist(replicate(K, S, simplify=FALSE)), dim=c(d, d, K))

  if (rotate) {
    res <- lol.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
  }
  # simulate from GMM
  sim <- lol.sims.sim_gmm(mus, S, n, priors)
  return(list(X=sim$X, Y=sim$Y, mus=mus, Sigmas=S))
}

#' Stacked Cigar
#'
#' A simulation for the stacked cigar experiment.
#' @import abind
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param rotate=FALSE whether to apply a random rotation to the mean and covariance. With random rotataion matrix Q, mu = Q*mu, and S = Q*S*Q.
#' @param priors=NULL the priors for each class. If NULL, class priors are all equal. If not null, should be |priors| = K, a length K vector for K classes.
#' @param a=0.15 scalar for all of the mu1 but 2nd dimension.
#' @param b=4 scalar for 2nd dimension value of mu2 and the 2nd variance term of S.
#' @return X [n, d] the data as a matrix.
#' @return Y [n] the labels as a array.
#' @return mus [d, K] the per-class means.
#' @return Sigmas [d, d, K] the per-class covariance matrices.
#' @author Eric Bridgeford
#' @examples
#' library(lol)
#' data <- lol.sims.cigar(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' @export
lol.sims.cigar <- function(n, d, rotate=FALSE, priors=NULL, a=0.15, b=4) {
  mu1 <- array(a, dim=c(d))
  mu1[2] <- b
  mus <- cbind(array(0, dim=c(d)), mu1)

  S <- diag(d)
  S[2,2] <- b

  S <- abind::abind(diag(d), S, along=3)

  if (rotate) {
    res <- lol.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
  }
  # simulate from GMM
  sim <- lol.sims.sim_gmm(mus, S, n, priors)
  return(list(X=sim$X, Y=sim$Y, mus=mus, Sigmas=S))
}


#' Xor Problem
#'
#' A function to simulate from the 2-class xor problem.
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param priors=NULL the priors for each class. If NULL, class priors are all equal. If not null, should be |priors| = K, a length K vector for K classes.
#' @param fall=100 the sigma for the covariance structuring. Sigma declines by ndim/fall across the variance terms.
#' @return X [n, d] the data as a matrix.
#' @return Y [n] the labels as a array.
#' @return mus [d, K] the per-class means.
#' @return Sigmas [d, d, K] the per-class covariance matrices.
#' @author Eric Bridgeford
#' @examples
#' library(lol)
#' data <- lol.sims.xor2(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' @export
lol.sims.xor2 <- function(n, d, priors=NULL, fall=100) {
  n1 <- ceiling(n/2)
  n2 <- floor(n/2)
  # first simulation set
  mus <- abind::abind(array(0, dim=c(d)), array(c(1, 0), dim=c(d)), along=2)
  S <- sqrt(d/fall)*diag(d)
  S <- abind::abind(S, S, along=3)

  # simulate from GMM for first set of training examples
  sim1 <- lol.sims.sim_gmm(mus, S, n1)

  # second simulation set
  mus <- abind::abind(array(1, dim=c(d)), array(c(0, 1), dim=c(d)), along=2)

  # simulate from GMM for second set of training examples
  sim2 <- lol.sims.sim_gmm(mus, S, n2)

  X <- abind::abind(sim1$X, sim2$X, along=1)
  Y <- abind::abind(sim1$Y, sim2$Y, along=1)

  reorder <- sample(n)
  return(list(X=X[reorder,], Y=Y[reorder], mus=mus, Sigmas=S))
}

#' GMM Simulate
#'
#' A helper function for simulating from Gaussian Mixture.
#' @param mus [d, K] the mus for each class.
#' @param Sigmas [d,d,K] the Sigmas for each class.
#' @param n the of examples.
#' @param priors=NULL the priors for each class. If NULL, class priors are all equal. If not null, should be |priors| = K, a length K vector for K classes.
#' @return X [n, d] the simulated data.
#' @return Y [n] the labels for each data point.
#' @author Eric Bridgeford
#' @import MASS
lol.sims.sim_gmm <- function(mus, Sigmas, n, priors=NULL) {
  K <- dim(mus)[2]
  if (is.null(priors)) {
    priors <- array(1/K, dim=c(K))
  } else if (length(priors) != K) {
    stop(sprintf("You have specified %d priors for %d classes.", length(priors), K))
  }
  labs <- sample(1:K, size=n, prob=priors, replace=TRUE)
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
lol.sims.rotation <- function(d) {
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
lol.sims.random_rotate <- function(mus, Sigmas, Q=NULL) {
  dimm <- dim(mus)
  K <- dimm[2]
  d <- dim(mus)[1]
  if (is.null(Q)) {
    Q <- lol.sims.rotation(d)
  } else if (!isTRUE(all.equal(dim(Q), c(d, d)))) {
    stop(sprintf("You have specified a rotation matrix with dimensions (%d, %d), but should be (%d, %d).", dim(Q)[1], dim(Q)[2], d, d))
  }

  for (i in 1:K) {
    mus[,i] <- Q %*% mus[,i,drop=FALSE]
    Sigmas[,,i] <- Q %*% Sigmas[,,i] %*% t(Q)
  }
  return(list(mus=mus, S=Sigmas, Q=Q))
}
