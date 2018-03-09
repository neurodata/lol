# Parallelize Stuff
#=========================#
require(MASS)
library(parallel)
require(lolR)
source('./plsda.R')

no_cores = detectCores() - 1

cl = makeCluster(no_cores)

# Setup Algorithms
#==========================#
algs <- list(lol.project.pca, lol.project.cpca, lol.project.lrcca, lol.project.pls, lol.project.lol, lol.project.qoq)
names(algs) <- c("PCA", "LDA", "CCA", "PLS", "LOL", "QOQ")
experiments <- list()
counter <- 1

# Setup Real Data
#==========================#
rlen <- 30
# the simulations to call themselves
maxr <- c(20, 30)
dat.names = c("Colon", "MNIST")

data <- list()

# Colon Data
path_colx <- 'data/colon/X.txt'
path_coly <- 'data/colon/Y.txt'
Xcol <- t(read.table(path_colx, sep="", header=FALSE, na.strings="", stringsAsFactors = FALSE))  # load in the data
Ycol <- as.vector(read.table(path_coly, sep="", header=FALSE, na.strings="", stringsAsFactors = FALSE))[,1]  # load in the labels as a 62-element vector
Ycol <- as.numeric(Ycol > 0)  # binarize for 2-class
data$Colon <- list(X=Xcol, Y=Ycol, exp="Colon")
rs <- unique(seq(from=1, to=maxr[1], length.out=min(maxr[1], rlen)))
for (r in rs) {
  for (j in 1:length(algs)) {
    alg <- algs[j]
    algname <- names(algs)[j]
    experiments[[counter]] <- list(exp="Colon", r=r, name=algname, k='loo', alg=alg)
    counter <- counter + 1
  }
}


# MNIST Data
load_mnist <- function(sourcedir='.') {
  load_image_file <- function(filename) {
    ret = list()
    f = file(filename,'rb')
    readBin(f,'integer',n=1,size=4,endian='big')
    ret$n = readBin(f,'integer',n=1,size=4,endian='big')
    nrow = readBin(f,'integer',n=1,size=4,endian='big')
    ncol = readBin(f,'integer',n=1,size=4,endian='big')
    x = readBin(f,'integer',n=ret$n*nrow*ncol,size=1,signed=F)
    ret$x = matrix(x, ncol=nrow*ncol, byrow=T)
    close(f)
    ret
  }
  load_label_file <- function(filename) {
    f = file(filename,'rb')
    readBin(f,'integer',n=1,size=4,endian='big')
    n = readBin(f,'integer',n=1,size=4,endian='big')
    y = readBin(f,'integer',n=n,size=1,signed=F)
    close(f)
    y
  }
  train <- list()
  test <- list()

  train$x <- load_image_file(paste(sourcedir, 'mnist/train-images-idx3-ubyte', sep='/'))
  test$x <- load_image_file(paste(sourcedir, 'mnist/t10k-images-idx3-ubyte', sep='/'))

  train$y <- load_label_file(paste(sourcedir, 'mnist/train-labels-idx1-ubyte', sep='/'))
  test$y <- load_label_file(paste(sourcedir, 'mnist/t10k-labels-idx1-ubyte', sep='/'))
  return(list(X.train=train$x, X.test=test$x, Y.train=train$y, Y.test=test$y))
}

show_digit <- function(png, title="",xlabel="Pixel", ylabel="Pixel", legend.name="metric", legend.show=TRUE,
                       font.size=12, limits=NULL) {
  mtx <- matrix(png, nrow=28)[,28:1]
  dm <- reshape2::melt(mtx)
  if (is.null(limits)) {
    limits <- c(min(mtx), max(mtx))
  }
  colnames(dm) <- c("x", "y", "value")
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  sqplot <- ggplot(dm, aes(x=x, y=y, fill=value)) +
    geom_tile() +
    scale_fill_gradientn(colours=jet.colors(7), name=legend.name, limits=limits) +
    xlab(xlabel) +
    ylab(ylabel) +
    ggtitle(title)
  if (legend.show) {
    sqplot <- sqplot +
      theme(text=element_text(size=font.size))
  } else {
    sqplot <- sqplot +
      theme(text=element_text(size=font.size, legend.position="none"))
  }
  return(sqplot)
}

result <- load_mnist('./data')
n <- 100  # the number of each class to train with
X <- array(0, dim=c(0, dim(result$X.train$x)[2]))
Y <- c()
for (lab in unique(result$Y.train)) {
  ss <- which(result$Y.train == lab)
  sn <- sample(ss, n)
  X <- rbind(X, result$X.train$x[sn,])
  Y <- c(Y, result$Y.train[sn])
}

res <- sample(1:length(Y), replace=FALSE)
X <- X[res,]; Y <- Y[res]  # randomly reorder the digits
data$MNIST <- list(X=X, Y=Y, exp="MNIST")
rs <- unique(seq(from=1, to=maxr[2], length.out=min(maxr[2], rlen)))
for (r in rs) {
  for (j in 1:length(algs)) {
    alg <- algs[j]
    algname <- names(algs)[j]
    experiments[[counter]] <- list(exp="MNIST", r=r, name=algname, k=10, alg=alg)
    counter <- counter + 1
  }
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
  }, error=function(e) lhat <- NaN)
  exr <- data.frame(data=exp$exp, se=sd(xv_res$Lhats)/sqrt(length(Y)), alg=exp$name, r=exp$r, lhat=lhat)
  return(exr)
})

# Aggregate and save
#=================================#
resultso <- do.call(rbind, results)
saveRDS(resultso, file.path(opath, 'lol_fig4_lda.rds'))
stopCluster(cl)
