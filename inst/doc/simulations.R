## ------------------------------------------------------------------------
require(lol)
require(ggplot2)
require(latex2exp)

n <- 1000
d <- 15
plot_sim <- function(X, Y, name) {
  data <- data.frame(x1=X[,1], x2=X[,2], y=Y)
  data$y <- factor(data$y)
  ggplot(data, aes(x=x1, y=x2, color=y)) +
    geom_point() +
    xlab("x1") +
    ylab("x2") +
    ggtitle(name)
}

## ---- fig.width=5--------------------------------------------------------
testdat <- lol.sims.rtrunk(n, d, b=20)
X <- testdat$X
Y <- testdat$Y
print(plot_sim(X, Y, "Trunk, 2 Class"))

## ---- fig.width=5--------------------------------------------------------
testdat <- lol.sims.rtrunk(n, d, rotate=TRUE, priors=c(0.8, 0.2), b=20)
X <- testdat$X
Y <- testdat$Y
print(plot_sim(X, Y, "Rotated Trunk, 2 Class, non-equal priors"))

## ---- fig.width=5--------------------------------------------------------
testdat <- lol.sims.rtrunk(n, d, b=20, K=3)
X <- testdat$X
Y <- testdat$Y
print(plot_sim(X, Y, "Trunk, 3 Class"))

## ---- fig.width=5--------------------------------------------------------
testdat <- lol.sims.mean_diff(n, d)
X <- testdat$X
Y <- testdat$Y
print(plot_sim(X, Y, "Mean Difference 2 Class"))

## ---- fig.width=5--------------------------------------------------------
testdat <- lol.sims.toep(n, d)
X <- testdat$X
Y <- testdat$Y
print(plot_sim(X, Y, "Toeplitz"))

## ---- fig.width=5--------------------------------------------------------
testdat <- lol.sims.qdtoep(n, d)
X <- testdat$X
Y <- testdat$Y
print(plot_sim(X, Y, "QD-Toeplitz"))

## ---- fig.width=5--------------------------------------------------------
testdat <- lol.sims.xor2(n, d)
X <- testdat$X
Y <- testdat$Y
print(plot_sim(X, Y, "XOR"))

## ---- fig.width=5--------------------------------------------------------
testdat <- lol.sims.cigar(n, d)
X <- testdat$X
Y <- testdat$Y
print(plot_sim(X, Y, "Cigar"))

## ---- fig.width=5--------------------------------------------------------
testdat <- lol.sims.fat_tails(n, d)
X <- testdat$X
Y <- testdat$Y
print(plot_sim(X, Y, "Fat Tails"))

