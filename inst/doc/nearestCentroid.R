## ---- echo=FALSE---------------------------------------------------------
require(lol)
require(ggplot2)
require(latex2exp)
require(MASS)
n=400
d=2

## ------------------------------------------------------------------------
testdat <- lol.sims.mean_diff(n, d)
X <- testdat$X
Y <- testdat$Y

data <- data.frame(x1=X[,1], x2=X[,2], y=Y)
data$y <- factor(data$y)
ggplot(data, aes(x=x1, y=x2, color=y)) +
  geom_point() +
  xlab("x1") +
  ylab("x2") +
  ggtitle("Simulated Data") +
  xlim(-4, 6) +
  ylim(-4, 4)

## ------------------------------------------------------------------------
classifier <- lol.classify.nearestCentroid(X, Y)

data <- cbind(data, data.frame(size=1))
data <- rbind(data, data.frame(x1=classifier$centroids[,1], x2=classifier$centroids[,2], y="center", size=5))
ggplot(data, aes(x=x1, y=x2, color=y, size=size)) +
  geom_point() +
  xlab("x1") +
  ylab("x2") +
  ggtitle("Data with estimated Centers") +
  guides(size=FALSE) +
  xlim(-4, 6) +
  ylim(-4, 4)

## ------------------------------------------------------------------------
Yhat <- predict(classifier, X)
data$y[1:(length(data$y) - 2)] <- Yhat
ggplot(data, aes(x=x1, y=x2, color=y, size=size)) +
  geom_point() +
  xlab("x1") +
  ylab("x2") +
  ggtitle("Data with Predictions") +
  guides(size=FALSE) +
  xlim(-4, 6) +
  ylim(-4, 4)

