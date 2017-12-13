#' Discriminant function which measures the distance between a point and a centroid j using Mahalanobis distance
#' @param x [1,r] data point that has dimensionality up to the rth subspace
#' @param centroid [1,r] centroid point j that has dimensionality up to the rth subspace
#' @param prior [1] prior probability for class j
#' @return Mahalanobis distance between x and centroid j
discriminant_fun <- function(x, centroid, prior){
  return(0.5 * sum((x - centroid)^2) - log(prior))
}

#' Linear Discriminant Analysis (LDA) Prediction
#'
#' A function for using Linear Discriminant Analysis (LDA) for prediction.
#' @param X [n, d] data matrix with n samples in d dimensions.
#' @param ylabs [C] vector containing the unique, ordered class labels.
#' @param centroids [C, d] centroid matrix of the unique, ordered class labels.
#' @param priors [C] vector containing prior probability for the unique, ordered class labels.
#' @param A [d, C-1] the projection matrix from d to C-1 dimensions.
#' @return Xr [n, C-1] projected data matrix.
#' @return Mp [C, C-1] projected centroid matrix.
#' @return Yhat [n] prediction matrix containing predictions for each example.
#' @author Richard Chen
#' @export
fs.predict.lda <- function(X, ylabs, centroids, priors, A){
  dimx <- dim(X)
  n <- dimx[1] # number of examples
  d <- dimx[2] # dimensionality of data
  C <- length(ylabs) # number of classes in the training set

  # Project the test data into the invariant subspaces
  Xr <- X %*% A
  # project the centroids into the invariant subspaces
  Mp <- centroids %*% A

  # Classify the data by doing nearest centroid classification
  Yhat <- ylabs[sapply(1:n, function(i) {
    which.min(sapply(1:C, function(j) {
      discriminant_fun(Xr[i,], Mp[j,], priors[j])
    }))
  })]

  return(list(Xr = Xr, Mp=Mp, Yhat=Yhat))
}

#' Low-rank Linear Discriminant Analysis (LR-LDA)
#'
#' A function for implementing the LR-LDA Algorithm.
#' @import Rspectra
#' @param X [n, d] the data with n samples in d dimensions.
#' @param Y [n] the labels of the samples.
#' @return A [d, C-1] the projection matrix
#' @return ylabs [C] vector containing the unique, ordered class labels.
#' @return centroids [C, d] centroid matrix of the unique, ordered classes.
#' @return priors [C] vector containing prior probability for the unique, ordered classes.
#' @author Richard Chen
#' @export
fs.project.lrlda <- function(X, Y) {
  classdat <- gs.utils.classdat(X, Y)
  priors <- classdat$priors; M <- classdat$centroids
  K <- classdat$C; ylabs <- classdat$ylabs
  n <- classdat$n; d <- classdat$d

  # 1. Compute the class dependent probabilities and class centroids
  priors = sapply(ylabs, function(y) sum(Y == y)/n)
  M <- t(as.matrix(sapply(ylabs, function(y) colMeans(X[Y==y,,drop=FALSE]))))

  # 2. Computing with-class covariance
  W <- cov(X) # within-class scatter

  # 3. Computing M* = M W^{-1/2} using the eigen-decomposition of W
  e <- eigen(W)
  V <- e$vectors # Recall that W decomposes to = V W V^T, which is V %*% diag(e$values) %*% t(V)
  W_neg_one_half <- V %*% diag(1./sqrt(e$values)) %*% t(V)
  M_star <- M %*% W_neg_one_half # M* = M W^{-1/2}

  # 4. Compuinge B* (which is just the covariance matrix of M*), and its eigen-decomposition, B*=V*D_BV*^T
  #    Note that the columns of V* define the coordiantes of the optimal subspaces
  B_star = cov(M_star)
  V_star = eigs(B_star, k=C-1)$vectors

  # 5. Compute the full projection matrix
  A = W_neg_one_half %*% V_star

  return(list(A=A, centroids=M, priors=priors, ylabs=ylabs))
}
