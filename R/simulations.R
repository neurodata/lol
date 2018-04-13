#' Fat Tails Simulation
#'
#' A function for simulating from 2 classes with differing means each with 2 sub-clusters, where one sub-cluster has a narrow tail and the other sub-cluster has a fat tail.
#'
#' @importFrom abind abind
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param rotate whether to apply a random rotation to the mean and covariance. With random rotataion matrix \code{Q}, \code{mu = Q*mu}, and \code{S = Q*S*Q}. Defaults to \code{FALSE}.
#' @param priors the priors for each class. If \code{NULL}, class priors are all equal. If not null, should be \code{|priors| = K}, a length \code{K} vector for \code{K} classes. Defaults to \code{NULL}.
#' @param f the fatness scaling of the tail. S2 = f*S1, where S1_{ij} = rho if i != j, and 1 if i == j. Defaults to \code{15}.
#' @param s0 the number of dimensions with a difference in the means. s0 should be < d. Defaults to \code{10}.
#' @param rho the scaling of the off-diagonal covariance terms, should be < 1. Defaults to \code{0.2}.
#' @param t the fraction of each class from the narrower-tailed distribution. Defaults to \code{0.8}.
#' @return A list of class \code{simulation} with the following:
#' \item{X}{\code{[n, d]} the \code{n} data points in \code{d} dimensions as a matrix.}
#' \item{Y}{\code{[n]} the \code{n} labels as an array.}
#' \item{mus}{\code{[d, K]} the \code{K} class means in \code{d} dimensions.}
#' \item{Sigmas}{\code{[d, d, K]} the \code{K} class covariance matrices in \code{d} dimensions.}
#' \item{priors}{\code{[K]} the priors for each of the \code{K} classes.}
#' \item{simtype}{The name of the simulation.}
#' \item{params}{Any extraneous parameters the simulation was created with.}
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("sims", package = "lolR")}
#'
#' @author Eric Bridgeford
#' @examples
#' library(lolR)
#' data <- lol.sims.fat_tails(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' @export
lol.sims.fat_tails <- function(n, d, rotate=FALSE, f=15, s0=10, rho=0.2, t=0.8, priors=NULL) {
  K <- 2
  if (is.null(priors)) {
    priors <- array(1/K, dim=c(K))
  } else if (length(priors) != K) {
    stop(sprintf("You have specified %d priors for %d classes.", length(priors), K))
  } else if (sum(priors) != 1) {
    stop(sprintf("You have passed invalid priors. The sum(priors) should be 1; yours is %.3f", sum(priors)))
  }
  if  (s0 > d) {
    stop(sprintf("s0 = %d, d=%d. s0 should be < d.", s0, d))
  }
  if (t > 1) {
    stop(sprintf("t = %.3f, while it should be a probability < 1.", t))
  }
  mu0 <- array(0, dim=c(d, 1))
  mu1 <- c(array(1, dim=c(s0)), array(0, dim=c(d - s0)))
  Q <- lol.sims.rotation(d)
  mus <- abind(mu0, mu1, along=2)

  S <- array(rho, dim=c(d, d))
  diag(S) <- 1

  S <- abind(S, S, along=3)

  if (rotate) {
    res <- lol.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
  }

  # simulate from GMM for narrow-tailed portion
  sim <- lol.sims.sim_gmm(mus, S, round(n*t), priors)
  # simulated from GMM for fat-tailed portion
  sim2 <- lol.sims.sim_gmm(mus, f*S, round(n*(1-t)), priors)
  X <- abind(sim$X, sim2$X, along=1)
  Y <- c(sim$Y, sim2$Y)
  return(structure(list(X=X, Y=Y, mus=mus, Sigmas=S, priors=sim$priors, simtype="Fat Tails",
                        params=list(f=f, s0=s0, rho=rho, t=t)), class="simulation"))
}

#' Mean Difference Simulation
#'
#' A function for simulating data in which a difference in the means is present only in a subset of dimensions, and equal covariance.
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param rotate whether to apply a random rotation to the mean and covariance. With random rotataion matrix \code{Q}, \code{mu = Q*mu}, and \code{S = Q*S*Q}. Defaults to \code{FALSE}.
#' @param priors the priors for each class. If \code{NULL}, class priors are all equal. If not null, should be \code{|priors| = K}, a length \code{K} vector for \code{K} classes. Defaults to \code{NULL}.
#' @param K the number of classes. Defaults to \code{2}.
#' @param md the magnitude of the difference in the means in the specified subset of dimensions. Ddefaults to \code{1}.
#' @param subset the dimensions to have a difference in the means. Defaults to only the first dimension. \code{max(subset) < d}. Defaults to \code{c(1)}.
#' @param offdiag the off-diagonal elements of the covariance matrix. Should be < 1. \code{S_{ij} = offdiag} if \code{i != j}, or 1 if \code{i == j}. Defaults to \code{0}.
#' @param s the scaling parameter of the covariance matrix. S_{ij} = scaling*1 if i == j, or scaling*offdiag if i != j. Defaults to \code{1}.
#' @return A list of class \code{simulation} with the following:
#' \item{X}{\code{[n, d]} the \code{n} data points in \code{d} dimensions as a matrix.}
#' \item{Y}{\code{[n]} the \code{n} labels as an array.}
#' \item{mus}{\code{[d, K]} the \code{K} class means in \code{d} dimensions.}
#' \item{Sigmas}{\code{[d, d, K]} the \code{K} class covariance matrices in \code{d} dimensions.}
#' \item{priors}{\code{[K]} the priors for each of the \code{K} classes.}
#' \item{simtype}{The name of the simulation.}
#' \item{params}{Any extraneous parameters the simulation was created with.}
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("sims", package = "lolR")}
#'
#' @author Eric Bridgeford
#' @examples
#' library(lolR)
#' data <- lol.sims.mean_diff(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' @export
lol.sims.mean_diff <- function(n, d, rotate=FALSE, priors=NULL, K=2, md=1, subset=c(1), offdiag=0, s=1) {
  if (is.null(priors)) {
    priors <- array(1/K, dim=c(K))
  } else if (length(priors) != K) {
    stop(sprintf("You have specified %d priors for %d classes.", length(priors), K))
  } else if (sum(priors) != 1) {
    stop(sprintf("You have passed invalid priors. The sum(priors) should be 1; yours is %.3f", sum(priors)))
  }
  if (max(subset) > d) {
    stop(sprintf("Specified a difference in dimension %d; maximum should be %d.", max(subset), d))
  }
  mus <- array(0, dim=c(d, K))
  for (i in 1:K) {
    mus[subset] <- (i - 1)*md  # scale the difference in the means accordingly
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
  return(structure(list(X=sim$X, Y=sim$Y, mus=mus, Sigmas=S, priors=sim$priors, simtype="Mean Difference",
                        params=list(K=K, md=md, subset=subset, offdiag=offdiag, s=s)), class="simulation"))
}

#' Toeplitz Simulation
#'
#' A function for simulating data in which the covariance is a non-symmetric toeplitz matrix.
#'
#' @importFrom abind abind
#' @importFrom stats toeplitz
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param rotate whether to apply a random rotation to the mean and covariance. With random rotataion matrix \code{Q}, \code{mu = Q*mu}, and \code{S = Q*S*Q}. Defaults to \code{FALSE}.
#' @param priors the priors for each class. If \code{NULL}, class priors are all equal. If not null, should be \code{|priors| = K}, a length \code{K} vector for \code{K} classes. Defaults to \code{NULL}.
#' @param D1 the dimensionality for the non-equal covariance terms. Defaults to \code{10}.
#' @param b a scaling parameter for the means. Defaults to \code{0.4}.
#' @param rho the scaling of the covariance terms, should be < 1. Defaults to \code{0.5}/
#' @return A list of class \code{simulation} with the following:
#' \item{X}{\code{[n, d]} the \code{n} data points in \code{d} dimensions as a matrix.}
#' \item{Y}{\code{[n]} the \code{n} labels as an array.}
#' \item{mus}{\code{[d, K]} the \code{K} class means in \code{d} dimensions.}
#' \item{Sigmas}{\code{[d, d, K]} the \code{K} class covariance matrices in \code{d} dimensions.}
#' \item{priors}{\code{[K]} the priors for each of the \code{K} classes.}
#' \item{simtype}{The name of the simulation.}
#' \item{params}{Any extraneous parameters the simulation was created with.}
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("sims", package = "lolR")}
#'
#' @author Eric Bridgeford
#' @examples
#' library(lolR)
#' data <- lol.sims.toep(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' @export
lol.sims.toep <- function(n, d, rotate=FALSE, priors=NULL, D1=10, b=0.4, rho=0.5) {
  K <- 2
  if (is.null(priors)) {
    priors <- array(1/K, dim=c(K))
  } else if (length(priors) != K) {
    stop(sprintf("You have specified %d priors for %d classes.", length(priors), K))
  } else if (sum(priors) != 1) {
    stop(sprintf("You have passed invalid priors. The sum(priors) should be 1; yours is %.3f", sum(priors)))
  }
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
  mus <- abind(mu0, -mu0, along=2)
  S <- abind(A, A, along=3)

  if (rotate) {
    res <- lol.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
  }
  # simulate from GMM
  sim <- lol.sims.sim_gmm(mus, S, n, priors)
  return(structure(list(X=sim$X, Y=sim$Y, mus=mus, Sigmas=S, priors=sim$priors, simtype="Toeplitz",
                        params=list(D1=D1, b=b, rho=rho)), class="simulation"))
}

#' Quadratic Discriminant Toeplitz Simulation
#'
#' A function for simulating data generalizing the Toeplitz setting, where each class has a different covariance matrix. This results in a Quadratic Discriminant.
#'
#' @importFrom abind abind
#' @importFrom stats toeplitz
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param rotate whether to apply a random rotation to the mean and covariance. With random rotataion matrix \code{Q}, \code{mu = Q*mu}, and \code{S = Q*S*Q}. Defaults to \code{FALSE}.
#' @param priors the priors for each class. If \code{NULL}, class priors are all equal. If not null, should be \code{|priors| = K}, a length \code{K} vector for \code{K} classes. Defaults to \code{NULL}.
#' @param D1 the dimensionality for the non-equal covariance terms. Defaults to \code{10}.
#' @param b a scaling parameter for the means. Defaults to \code{0.4}.
#' @param rho the scaling of the covariance terms, should be < 1. Defaults to \code{0.5}.
#' @return A list of class \code{simulation} with the following:
#' \item{X}{\code{[n, d]} the \code{n} data points in \code{d} dimensions as a matrix.}
#' \item{Y}{\code{[n]} the \code{n} labels as an array.}
#' \item{mus}{\code{[d, K]} the \code{K} class means in \code{d} dimensions.}
#' \item{Sigmas}{\code{[d, d, K]} the \code{K} class covariance matrices in \code{d} dimensions.}
#' \item{priors}{\code{[K]} the priors for each of the \code{K} classes.}
#' \item{simtype}{The name of the simulation.}
#' \item{params}{Any extraneous parameters the simulation was created with.}
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("sims", package = "lolR")}
#'
#' @author Eric Bridgeford
#' @examples
#' library(lolR)
#' data <- lol.sims.qdtoep(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' @export
lol.sims.qdtoep <- function(n, d, rotate=FALSE, priors=NULL, D1=10, b=0.4, rho=0.5) {
  K <- 2
  if (is.null(priors)) {
    priors <- array(1/K, dim=c(K))
  } else if (length(priors) != K) {
    stop(sprintf("You have specified %d priors for %d classes.", length(priors), K))
  } else if (sum(priors) != 1) {
    stop(sprintf("You have passed invalid priors. The sum(priors) should be 1; yours is %.3f", sum(priors)))
  }
  if (rho >= 1) {
    stop(sprintf("rho should be < 1; user specified %.3f", rho))
  }
  if (D1 > d) {
    stop(sprintf("User has specified more dimensions for non-equal cov terms. D1 is %d, yet d is %d", D1, d))
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
  mus <- abind(mu0, mu1, along=2)
  A1 <- Q %*% A0 %*% t(Q)
  S <- abind(A0, A1, along=3)

  if (rotate) {
    res <- lol.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
  }
  # simulate from GMM
  sim <- lol.sims.sim_gmm(mus, S, n, priors)
  return(structure(list(X=sim$X, Y=sim$Y, mus=mus, Sigmas=S, priors=sim$priors, simtype="Quadratic Toeplitz",
                        params=list(D1=D1, b=b, rho=rho)), class="simulation"))
}

#' Random Trunk
#'
#' A simulation for the random trunk experiment.
#' @importFrom abind abind
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param rotate whether to apply a random rotation to the mean and covariance. With random rotataion matrix \code{Q}, \code{mu = Q*mu}, and \code{S = Q*S*Q}. Defaults to \code{FALSE}.
#' @param priors the priors for each class. If \code{NULL}, class priors are all equal. If not null, should be \code{|priors| = K}, a length \code{K} vector for \code{K} classes. Defaults to \code{NULL}.
#' @param b scalar for mu scaling. Default to \code{4}.
#' @param K number of classes, should be <4. Defaults to \code{2}.
#' @return A list of class \code{simulation} with the following:
#' \item{X}{\code{[n, d]} the \code{n} data points in \code{d} dimensions as a matrix.}
#' \item{Y}{\code{[n]} the \code{n} labels as an array.}
#' \item{mus}{\code{[d, K]} the \code{K} class means in \code{d} dimensions.}
#' \item{Sigmas}{\code{[d, d, K]} the \code{K} class covariance matrices in \code{d} dimensions.}
#' \item{priors}{\code{[K]} the priors for each of the \code{K} classes.}
#' \item{simtype}{The name of the simulation.}
#' \item{params}{Any extraneous parameters the simulation was created with.}
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("sims", package = "lolR")}
#'
#' @author Eric Bridgeford
#' @examples
#' library(lolR)
#' data <- lol.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' @export
lol.sims.rtrunk <- function(n, d, rotate=FALSE, priors=NULL, b=4, K=2) {
  if (is.null(priors)) {
    priors <- array(1/K, dim=c(K))
  } else if (length(priors) != K) {
    stop(sprintf("You have specified %d priors for %d classes.", length(priors), K))
  } else if (sum(priors) != 1) {
    stop(sprintf("You have passed invalid priors. The sum(priors) should be 1; yours is %.3f", sum(priors)))
  }
  mu1 <- b/sqrt(0:(d-1)*2 + 1)
  if (K == 2) {
    mus <- abind(mu1, -mu1, along=2)
  } else if (K == 3) {
    mus <- abind(mu1, 0*mu1, -mu1, along=2)
  } else if (K == 4) {
    mus <- abind(mu1, b/(seq(from=d, to=1, by=-1)), b/(seq(from=1, to=d, by=1)), -mu1, along=2)
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
  return(structure(list(X=sim$X, Y=sim$Y, mus=mus, Sigmas=S, priors=sim$priors, simtype="Random Trunk",
                        params=list(b=b, K=K)), class="simulation"))
}

#' Stacked Cigar
#'
#' A simulation for the stacked cigar experiment.
#' @importFrom abind abind
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param rotate whether to apply a random rotation to the mean and covariance. With random rotataion matrix \code{Q}, \code{mu = Q*mu}, and \code{S = Q*S*Q}. Defaults to \code{FALSE}.
#' @param priors the priors for each class. If \code{NULL}, class priors are all equal. If not null, should be \code{|priors| = K}, a length \code{K} vector for \code{K} classes. Defaults to \code{NULL}.
#' @param a scalar for all of the mu1 but 2nd dimension. Defaults to \code{0.15}.
#' @param b scalar for 2nd dimension value of mu2 and the 2nd variance term of S. Defaults to \code{4}.
#' @return A list of class \code{simulation} with the following:
#' \item{X}{\code{[n, d]} the \code{n} data points in \code{d} dimensions as a matrix.}
#' \item{Y}{\code{[n]} the \code{n} labels as an array.}
#' \item{mus}{\code{[d, K]} the \code{K} class means in \code{d} dimensions.}
#' \item{Sigmas}{\code{[d, d, K]} the \code{K} class covariance matrices in \code{d} dimensions.}
#' \item{priors}{\code{[K]} the priors for each of the \code{K} classes.}
#' \item{simtype}{The name of the simulation.}
#' \item{params}{Any extraneous parameters the simulation was created with.}
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("sims", package = "lolR")}
#'
#' @author Eric Bridgeford
#' @examples
#' library(lolR)
#' data <- lol.sims.cigar(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' @export
lol.sims.cigar <- function(n, d, rotate=FALSE, priors=NULL, a=0.15, b=4) {
  K <- 2
  if (is.null(priors)) {
    priors <- array(1/K, dim=c(K))
  } else if (length(priors) != K) {
    stop(sprintf("You have specified %d priors for %d classes.", length(priors), K))
  } else if (sum(priors) != 1) {
    stop(sprintf("You have passed invalid priors. The sum(priors) should be 1; yours is %.3f", sum(priors)))
  }
  mu1 <- array(a, dim=c(d))
  mu1[2] <- b
  mus <- cbind(array(0, dim=c(d)), mu1)

  S <- diag(d)
  S[2,2] <- b

  S <- abind(diag(d), S, along=3)

  if (rotate) {
    res <- lol.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
  }
  # simulate from GMM
  sim <- lol.sims.sim_gmm(mus, S, n, priors)
  return(structure(list(X=sim$X, Y=sim$Y, mus=mus, Sigmas=S, priors=sim$priors, simtype="Stacked Cigar",
                        params=c(a=a, b=b)), class="simulation"))
}


#' Xor Problem
#'
#' A function to simulate from the 2-class xor problem.
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param priors the priors for each class. If \code{NULL}, class priors are all equal. If not null, should be \code{|priors| = K}, a length \code{K} vector for \code{K} classes. Defaults to \code{NULL}.
#' @param fall the falloff for the covariance structuring. Sigma declines by ndim/fall across the variance terms. Defaults to \code{100}.
#' @return A list of class \code{simulation} with the following:
#' \item{X}{\code{[n, d]} the \code{n} data points in \code{d} dimensions as a matrix.}
#' \item{Y}{\code{[n]} the \code{n} labels as an array.}
#' \item{mus}{\code{[d, K]} the \code{K} class means in \code{d} dimensions.}
#' \item{Sigmas}{\code{[d, d, K]} the \code{K} class covariance matrices in \code{d} dimensions.}
#' \item{priors}{\code{[K]} the priors for each of the \code{K} classes.}
#' \item{simtype}{The name of the simulation.}
#' \item{params}{Any extraneous parameters the simulation was created with.}
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("sims", package = "lolR")}
#'
#' @author Eric Bridgeford
#' @examples
#' library(lolR)
#' data <- lol.sims.xor2(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' @export
lol.sims.xor2 <- function(n, d, priors=NULL, fall=100) {
  K <- 2
  if (is.null(priors)) {
    priors <- array(1/K, dim=c(K))
  } else if (length(priors) != K) {
    stop(sprintf("You have specified %d priors for %d classes.", length(priors), K))
  } else if (sum(priors) != 1) {
    stop(sprintf("You have passed invalid priors. The sum(priors) should be 1; yours is %.3f", sum(priors)))
  }
  n1 <- ceiling(n/2)
  n2 <- floor(n/2)
  # first simulation set
  mus <- abind(array(0, dim=c(d)), array(c(1, 0), dim=c(d)), along=2)
  S <- sqrt(d/fall)*diag(d)
  S <- abind(S, S, along=3)

  # simulate from GMM for first set of training examples
  sim1 <- lol.sims.sim_gmm(mus, S, n1, priors)

  # second simulation set
  mus <- abind(array(1, dim=c(d)), array(c(0, 1), dim=c(d)), along=2)

  # simulate from GMM for second set of training examples
  sim2 <- lol.sims.sim_gmm(mus, S, n2, priors=priors)

  X <- abind(sim1$X, sim2$X, along=1)
  Y <- abind(sim1$Y, sim2$Y, along=1)

  reorder <- sample(n)
  return(structure(list(X=X[reorder,], Y=Y[reorder], mus=mus, Sigmas=S, priors=sim2$priors, simtype="Xor",
                        params=list(fall=fall)), class="simulation"))
}

#' GMM Simulate
#'
#' A helper function for simulating from Gaussian Mixture.
#' @param mus \code{[d, K]} the mus for each class.
#' @param Sigmas \code{[d,d,K]} the Sigmas for each class.
#' @param n the number of examples.
#' @param priors \code{K} the priors for each class.
#' @return A list with the following:
#' \item{X}{\code{[n, d]} the simulated data.}
#' \item{Y}{\code{[n]} the labels for each data point.}
#' \item{priors}{\code{[K]} the priors for each class.}
#' @author Eric Bridgeford
#' @import MASS
lol.sims.sim_gmm <- function(mus, Sigmas, n, priors) {
  K <- dim(mus)[2]
  labs <- sample(1:K, size=n, prob=priors, replace=TRUE)
  ylabs <- as.vector(sort(unique(labs)))
  res <- sapply(ylabs, function(y) mvrnorm(n=sum(labs == y), mus[,y], Sigmas[,,y]), USE.NAMES=TRUE, simplify=FALSE)
  X <- array(0, dim=c(n, dim(Sigmas)[1]))
  for (y in ylabs) {
    X[labs == y,] <- res[[y]]
  }
  return(list(X=X, Y=labs, priors=priors))
}

#' Sample Random Rotation
#'
#' A helper function for estimating a random rotation matrix.
#' @importFrom stats rnorm
#' @param d dimensions to generate a rotation matrix for.
#' @return the rotation matrix
#' @author Eric Bridgeford
lol.sims.rotation <- function(d) {
  Q <- qr.Q(qr(array(rnorm(d*d), dim=c(d, d))))
  if (det(Q) < -.99) {
    Q[,1] <- -Q[,1]
  }
  return(Q)
}

#' Random Rotation
#'
#' A helper function for applying a random rotation to gaussian parameter set.
#' @param mus means per class.
#' @param Sigmas covariances per class.
#' @param Q rotation to use, if any
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
