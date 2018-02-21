# Parallelize Stuff
#=========================#
require(MASS)
library(parallel)
require(lolR)
source('../../../datasets/PMLB/load.R')
source('../../../datasets/prepare_dataset.R')

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
pmlb.path <- '../../../penn-ml-benchmarks/'
dset.file <- file.path(pmlb.path, 'datasets', 'classification', 'classification_datasets_pmlb.tsv')
dset.names <- read.table(dset.file, sep="", header=TRUE)$name
rlen <- 30
cutoff <- 1000

data <- list()

log.seq <- function(from=0, to=30, length=15) {
  round(exp(seq(from=log(from), to=log(to), length.out=length)))
}
for (i in 1:length(dset.names)) {
  tryCatch({result <- load.pmlb(dset.names[i], pmlb.path, 'classification')
  result <- clean.dataset(result$X, result$Y, Kmax=10)
  data[[dset.names[i]]] <- list(X=result$X, Y=result$Y, exp=dset.names[i])
  n <- length(result$Y)
  if (n > cutoff) {
    k <- 10
  } else {
    k <- 'loo'
  }
  n <- dim(result$X)[1]; d <- dim(result$X)[2]
  maxr <- min(d, 100)
  rs <- unique(log.seq(from=1, to=maxr, length=15))
  for (r in rs) {
    for (j in 1:length(algs)) {
      alg <- algs[j]
      algname <- names(algs)[j]
      experiments[[counter]] <- list(exp=dset.names[i], r=r, name=algname, k=k, alg=alg)
      counter <- counter + 1
    }
  }}, error = function(e) NaN)
}

# Setup Algorithms
#=========================#
clusterExport(cl, "data"); clusterExport(cl, "rlen")
clusterExport(cl, "experiments")
results <- parLapply(cl, experiments, function(exp) {
  require(lolR)
  alg <- exp$alg
  name <- exp$name
  results <- data.frame(exp=c(), alg=c(), r=c(), lhat=c())
  if (name %in% c("QOQ")) {
    classifier.alg=MASS::qda
  } else {
    classifier.alg=MASS::lda
  }
  X <- data[[exp$exp]]$X; Y <- data[[exp$exp]]$Y
  tryCatch({
    xv_res <- lol.xval.eval(X, Y, alg=exp$alg[[exp$name]], alg.opts=list(r=exp$r), alg.embedding="A",
                            classifier=classifier.alg, k=exp$k)
    lhat <- xv_res$Lhat
    exr <- data.frame(data=exp$exp, se=sd(xv_res$Lhats)/sqrt(length(Y)), alg=exp$name, r=exp$r, K=length(unique(Y)), n = dim(X)[1], lhat=lhat)
  }, error=function(e) lhat <- NaN)
  return(exr)
})

# Aggregate and save
#=================================#
resultso <- do.call(rbind, results)
saveRDS(resultso, 'lol_fig4pmlb_lda.rds')
stopCluster(cl)
