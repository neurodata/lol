require(lolR)
require(FlashR)
require(abind)
source('../../flashlol-figure/flashLol.R')
## Simulation Tests

# All work on vertically stacked cigars
sim <- lol.sims.cigar(n=100, d=5)
X.cig <- fm.as.matrix(sim$X); Y.cig <- sim$Y
cig.pca <- flashx.pca(X=X.cig, Y=Y.cig, r=2)
cig.lol <- flashx.lol(X=X.cig, Y=Y.cig, r=2)
cig.lrlda <- flashx.lrlda(X=X.cig, Y=Y.cig, r=2)
cig.cca <- flashx.lrcca(X=X.cig, Y=Y.cig, r=2)
cig.rp <- flashx.rp(X=X.cig, Y=Y.cig, r=2)
cig.out <- list(sim=list(X=as.matrix(X.cig), Y=Y.cig), pca=cig.pca, lol=cig.lol,
                 lrlda=cig.lrlda, cca=cig.cca, rp=cig.rp)

# Horizontally stacked cigars breaks PCA
lol.sims.hcigar <- function(n, d, rotate=FALSE, priors=NULL, a=0.15, b=4) {
  K <- 2
  if (is.null(priors)) {
    priors <- array(1/K, dim=c(K))
  } else if (length(priors) != K) {
    stop(sprintf("You have specified %d priors for %d classes.", length(priors), K))
  } else if (sum(priors) != 1) {
    stop(sprintf("You have passed invalid priors. The sum(priors) should be 1; yours is %.3f", sum(priors)))
  }
  mu1 <- array(a, dim=c(d))
  mu1[1] <- b
  mus <- cbind(array(0, dim=c(d)), mu1)

  S <- diag(d)
  S[2,2] <- b^2

  S <- abind(S, S, along=3)

  if (rotate) {
    res <- lolR:::lol.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
  }
  # simulate from GMM
  sim <- lolR:::lol.sims.sim_gmm(mus, S, n, priors)
  return(structure(list(X=sim$X, Y=sim$Y, mus=mus, Sigmas=S, priors=sim$priors, simtype="Stacked Cigar",
                        params=c(a=a, b=b)), class="simulation"))
}
sim <- lol.sims.hcigar(n=100, d=5)
X.hcig <- fm.as.matrix(sim$X); Y.hcig <- sim$Y
hcig.pca <- flashx.pca(X=X.hcig, Y=Y.hcig, r=2)
hcig.lol <- flashx.lol(X=X.hcig, Y=Y.hcig, r=2)
hcig.lrlda <- flashx.lrlda(X=X.hcig, Y=Y.hcig, r=2)
hcig.cca <- flashx.lrcca(X=X.hcig, Y=Y.hcig, r=2)
hcig.rp <- flashx.rp(X=X.hcig, Y=Y.hcig, r=2)
hcig.out <- list(sim=list(X=as.matrix(X.hcig), Y=Y.hcig), pca=hcig.pca, lol=hcig.lol,
                 lrlda=hcig.lrlda, cca=hcig.cca, rp=hcig.rp)

# CCA Maximally Piles when d > n
sim <- lol.sims.rtrunk(n=100, d=101)
X.rt <- fm.as.matrix(sim$X); Y.rt <- sim$Y
rt.pca <- flashx.pca(X=X.rt, Y=Y.rt, r=2)
rt.lol <- flashx.lol(X=X.rt, Y=Y.rt, r=2)
rt.lrlda <- flashx.lrlda(X=X.rt, Y=Y.rt, r=2)
rt.cca <- flashx.lrcca(X=X.rt, Y=Y.rt, r=2)
rt.rp <- flashx.rp(X=X.rt, Y=Y.rt, r=2)
rt.out <- list(sim=list(X=as.matrix(X.rt), Y=Y.rt), pca=rt.pca, lol=rt.lol,
                 lrlda=rt.lrlda, cca=rt.cca, rp=rt.rp)

saveRDS(list(cig=cig.out, hcig=hcig.out, lowrank=rt.out), './flashx_sims.rds')
