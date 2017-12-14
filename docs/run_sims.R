nsim <- 10

ds <- c(15, 50, 100, 500)
n <- 100
nd <- length(ds)

r <- 10
files = c('xor_sim.rds', 'cig_sim.rds', 'rtr_sim.rds', 'rtr3_sim.rds', 'toe_sim.rds', 'ft_sim.rds')

performance <- data.frame(dimensions=c(), algorithm=c(), classification=c(), simulation=c(), simid=c(), lhat=c())
sims = c('xor', 'cigar', 'rtrunk', 'rtrunk3', 'toep', 'fat tails')
algorithms = c(fs.project.cpca, fs.project.lol, fs.project.pca, fs.project.lrcca)
algnames =c("cpca", "lol", "pca", "lrcca")
classalgs <- c("lda", "rf")

for (i in 3:length(sims)) {
  print(paste('Simulations:', sims[i]))
  simset <- readRDS(files[i])
  Xsets <- simset$X
  Ysets <- simset$Y
  for (j in 1:length(ds)) {
    print(paste('Dimensions:', ds[j]))
    Xd <- Xsets[[j]]
    Yd <- Ysets[[j]]
    for (k in 1:nsim) {
      cat(paste(k))
      X <- Xd[k,,]
      Y <- Yd[k,]
      for (l in 1:length(algorithms)) {
        for (m in 1:length(classalgs)) {
          res <- suppressWarnings(fs.eval.xval(X, Y, r=r, alg=algorithms[l][[1]], classifier=classalgs[m], k='loo'))
          performance <- rbind(performance, data.frame(dimensions=ds[j], algorithm=algnames[l], classification=classalgs[m],
                                                       simulation=sims[i], simid=k, lhat=res$Lhat))
        }
      }
      cat("\n")
    }
  }
}

saveRDS(performance, 'performance.rds')

