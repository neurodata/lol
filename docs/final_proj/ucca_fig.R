p <- 300
q <- 3
delta <- 0.4
S <- diag(p)
Sigmas <- abind(S, S, S, along=3)
mus <- array(NaN, dim=c(p, q))
for (j in 1:q) {
  mus[,j] <- delta*(j - 1)
}

sample <- lol:::lol.sims.sim_gmm(mus, Sigmas, 75, priors=rep(1/q, q))
test <- lol:::lol.sims.sim_gmm(mus, Sigmas, 600, priors=rep(1/q, q))

res <- lol.project.lrcca(sample$X, sample$Y, 2)
t.Xr <- test$X %*% res$A

liney <- lda(res$Xr, sample$Y)
result <- predict(liney, t.Xr)


data <- data.frame(x1=t.Xr[,1], x2=t.Xr[,2], y=test$Y)
data$y <- factor(data$y)

# Same, but with different shapes
ggplot(data, aes(x=x1, y=x2, shape=y, color=y)) + geom_point() +
  scale_shape_manual(values=c(1,2,5)) +
  ylim(-3, 3) +
  xlab(TeX("$X^Ta_1$")) +
  ylab(TeX("$X^Ta_2$"))

data <- data.frame(x1=result$x[,1], y=sample$Y)
data$y <- factor(data$y)
ggplot(data, aes(x=x1, color=y, group=y)) +
  geom_density()



data <- data.frame(x1=result$x[,1], x2=result$x[,2], y=sample$Y)
data$y <- factor(data$y)
ggplot(data, aes(x=x1, y=x2, color=y, group=y)) +
  geom_point()


data <- data.frame(x1=sample$X[,1], x2=sample$X[,2], y=sample$Y)
data$y <- factor(data$y)
ggplot(data, aes(x=x1, y=x2, color=y, group=y)) +
  geom_point()




data <- data.frame(x1=t.Xr[,1], x2=t.Xr[,2], y=sample$Y)
data$y <- factor(data$y)
ggplot(data, aes(x=x1, y=x2, color=y, group=y)) +
  geom_point()
