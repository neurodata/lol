# Parallelize Stuff
#=========================#
require(MASS)
library(parallel)
require(lolR)
require(slbR)

no_cores = detectCores() - 1

cl = makeCluster(no_cores)

# Setup Algorithms
#==========================#
algs <- list(lol.project.pca, lol.project.cpca, lol.project.lrcca, lol.project.lol, lol.project.qoq)
names(algs) <- c("PCA", "LDA", "CCA", "LOL", "QOQ")
experiments <- list()
counter <- 1

# Setup Real Data
#==========================#
rlen <- 15
ncutoff <- 1000

data <- list()

log.seq <- function(from=0, to=30, length=15) {
  round(exp(seq(from=log(from), to=log(to), length.out=length)))
}

dset.names <- names(pmlb.list(task="classification")$dsets.info)
for (i in 1:length(dset.names)) {
  tryCatch({result <- pmlb.load(datasets = dset.names[i], tasks='classification', clean.nan=TRUE, clean.ohe=TRUE)
  result <- result$data[[dset.names[i]]]
  data[[dset.names[i]]] <- list(X=result$X, Y=result$Y, exp=dset.names[i])
  n <- length(result$Y)
  if (n > ncutoff) {
    k <- 10
  } else {
    k <- 'loo'
  }
  for (j in 1:length(algs)) {
    alg <- algs[j]
    algname <- names(algs)[j]
    experiments[[counter]] <- list(exp=dset.names[i], name=algname, k=k, alg=alg)
    counter <- counter + 1
  }}, error = function(e) NaN)
}

# Setup Algorithms
#=========================#
opath <- '/home/ebridge/research/R-fmri/lol/docs/lol-paper/data/fig5/'
clusterExport(cl, "data"); clusterExport(cl, "rlen")
clusterExport(cl, "experiments"); clusterExport(cl, "opath")
results <- parLapply(cl, experiments, function(exp) {
  require(lolR)
  if (exp$name %in% c("QOQ")) {
    classifier.alg=MASS::qda
  } else {
    classifier.alg=MASS::lda
  }
  X <- data[[exp$exp]]$X; Y <- data[[exp$exp]]$Y
  n <- dim(X)[1]; d <- dim(X)[2]
  maxr <- min(d, 100)
  rs <- unique(log.seq(from=1, to=maxr, length=rlen))
  results <- data.frame(exp=c(), alg=c(), K=c(), r=c(), n=c(), lhat=c())
  for (r in rs) {
    tryCatch({
      xv_res <- lol.xval.eval(X, Y, alg=exp$alg[[exp$name]], alg.opts=list(r=r), alg.embedding="A",
                              classifier=classifier.alg, k=exp$k)
      lhat <- xv_res$Lhat
      exr <- data.frame(data=exp$exp, se=sd(xv_res$Lhats)/sqrt(length(Y)), alg=exp$name, r=r,
                        K=length(unique(Y)), n=n, lhat=lhat)
      results <- rbind(results, exr)
    }, error=function(e) lhat <- NaN)
  }
  saveRDS(results, file=paste(opath, exp$exp, '.rds', sep=""))
})
saveRDS(resultso, 'lol_fig4pmlb_lda.rds')
stopCluster(cl)

# Aggregate and save
#=================================#
results <- lapply(dset.names, function(dset) {
  result <- readRDS(paste(opath, dset, '.rds', sep=""))
})

resultso <- do.call(rbind, results)
