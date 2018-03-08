## ------------------------------------------------------------------------
require(lol)
require(ggplot2)
require(latex2exp)
n = 400
d = 30
r = 3

## ---- fig.width=5--------------------------------------------------------
testdat <- lol.sims.rtrunk(n, d)
X <- testdat$X
Y <- testdat$Y

data <- data.frame(x1=X[,1], x2=X[,2], y=Y)
data$y <- factor(data$y)
ggplot(data, aes(x=x1, y=x2, color=y)) +
  geom_point() +
  xlab("x1") +
  ylab("x2") +
  ggtitle("Simulated Data")

## ---- fig.width=5--------------------------------------------------------
result <- lol.xval.eval(X, Y, alg = lol.project.lol, alg.opts=list(r=r), alg.return="A",
                        classifier=MASS::lda, classifier.return="class", k='loo')

data <- data.frame(x1=result$model$Xr[,1], x2=result$model$Xr[,2], y=Y)
data$y <- factor(data$y)
ggplot(data, aes(x=x1, y=x2, color=y)) +
  geom_point() +
  xlab("x1") +
  ylab("x2") +
  ggtitle(sprintf("Projected Data using LOL, L=%.2f", result$Lhat))

