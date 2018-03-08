## ---- echo=FALSE---------------------------------------------------------
require(lol)
require(ggplot2)
require(latex2exp)
require(MASS)
n=400
d=30
r=3

## ---- fig.width=5--------------------------------------------------------
testdat <- lol.sims.cigar(n, d)
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
result <- lol.project.cpca(X, Y, r)

data <- data.frame(x1=result$Xr[,1], x2=result$Xr[,2], y=Y)
data$y <- factor(data$y)
ggplot(data, aes(x=x1, y=x2, color=y)) +
  geom_point() +
  xlab("x1") +
  ylab("x2") +
  ggtitle("Projected Data using cPCA")

## ---- fig.width=5--------------------------------------------------------
liney <- MASS::lda(result$Xr, Y)
result <- predict(liney, result$Xr)
lhat <- 1 - sum(result$class == Y)/length(Y)

data <- data.frame(x1=result$x[,1], y=Y)
data$y <- factor(data$y)
ggplot(data, aes(x=x1, fill=y)) +
  geom_density(adjust=1.5, alpha=0.6) +
  xlab("x1") +
  ylab("Density") +
  ggtitle(sprintf("cPCA-LDA, L = %.2f", lhat))

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
result <- lol.project.cpca(X, Y, r)

data <- data.frame(x1=result$Xr[,1], x2=result$Xr[,2], y=Y)
data$y <- factor(data$y)
ggplot(data, aes(x=x1, y=x2, color=y)) +
  geom_point() +
  xlab("x1") +
  ylab("x2") +
  ggtitle("Projected Data using cPCA")

## ---- fig.width=5--------------------------------------------------------
liney <- MASS::lda(result$Xr, Y)
result <- predict(liney, result$Xr)
lhat <- 1 - sum(result$class == Y)/length(Y)

data <- data.frame(x1=result$x[,1], y=Y)
data$y <- factor(data$y)
ggplot(data, aes(x=x1, fill=y)) +
  geom_density(adjust=1.5, alpha=0.6) +
  xlab("x1") +
  ylab("Density") +
  ggtitle(sprintf("cPA-LDA, L = %.2f", lhat))

## ---- fig.width=5--------------------------------------------------------
testdat <- lol.sims.rtrunk(n, d, rotate=TRUE)
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
result <- lol.project.cpca(X, Y, r)

data <- data.frame(x1=result$Xr[,1], x2=result$Xr[,2], y=Y)
data$y <- factor(data$y)
ggplot(data, aes(x=x1, y=x2, color=y)) +
  geom_point() +
  xlab("x1") +
  ylab("x2") +
  ggtitle("Projected Data using cPCA")

## ---- fig.width=5--------------------------------------------------------
liney <- MASS::lda(result$Xr, Y)
result <- predict(liney, result$Xr)
lhat <- 1 - sum(result$class == Y)/length(Y)

data <- data.frame(x1=result$x[,1], y=Y)
data$y <- factor(data$y)
ggplot(data, aes(x=x1, fill=y)) +
  geom_density(adjust=1.5, alpha=0.6) +
  xlab("x1") +
  ylab("Density") +
  ggtitle(sprintf("cPCA-LDA, L = %.2f", lhat))

