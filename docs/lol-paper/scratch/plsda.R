require(caret)
require(pls)
require(R.utils)

lol.project.pls <- function(X, Y, r, ...) {
  info <- lol:::lol.utils.info(X, Y)
  priors <- info$priors; centroids <- info$centroids
  K <- info$K; ylabs <- info$ylabs
  n <- info$n; d <- info$d
  A <- plsda(X, as.factor(Y), r)$projection
  return(list(A=A, Xa=lol.embed(X, A), ylabs=ylabs))
}
