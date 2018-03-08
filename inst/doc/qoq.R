## ---- echo=FALSE---------------------------------------------------------
require(lol)
require(ggplot2)
require(latex2exp)
require(MASS)
n=100
d=100
r=5

## ---- fig.width=5--------------------------------------------------------
testdat <- lol.sims.qdtoep(n, d)
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
result <- lol.project.qoq(X, Y, r)

data <- data.frame(x1=result$Xr[,1], x2=result$Xr[,2], y=Y)
data$y <- factor(data$y)
ggplot(data, aes(x=x1, y=x2, color=y)) +
  geom_point() +
  xlab("x1$") +
  ylab("x2") +
  ggtitle("Projected Data using QOQ")

## ---- fig.width=5--------------------------------------------------------
quaddy <- MASS::qda(result$Xr, Y)
result <- predict(quaddy, result$Xr)
lhat <- 1 - sum(result$class == Y)/length(Y)

print(sprintf("QOQ, QDA L =%.3f", lhat))

## ------------------------------------------------------------------------
resultl <- lol.project.lol(X, Y, r)

data <- data.frame(x1=resultl$Xr[,1], x2=resultl$Xr[,2], y=Y)
data$y <- factor(data$y)
ggplot(data, aes(x=x1, y=x2, color=y)) +
  geom_point() +
  xlab("x1") +
  ylab("x2") +
  ggtitle("Projected Data using LOL")

## ------------------------------------------------------------------------
liney <- MASS::qda(resultl$Xr, Y)
result <- predict(liney, resultl$Xr)
lhat <- 1 - sum(result$class == Y)/length(Y)

print(sprintf("LOL, LDA L =%.3f", lhat))

quaddy <- MASS::qda(resultl$Xr, Y)
result <- predict(quaddy, resultl$Xr)
lhat <- 1 - sum(result$class == Y)/length(Y)

print(sprintf("LOL, QDA  L =%.3f", lhat))

