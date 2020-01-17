## Projection Analysis
require(oro.nifti)
require(ramify)
require(parallel)
require(profmem)
require(FlashR)
require(MASS)
require(dplyr)
source('flashLol.R')

mc.cores <- detectCores()
fnames <- list.files('/brains/brains')
fnames <- fnames[fnames != "lost+found"]

subids <- as.integer(sapply(strsplit(fnames, "-|_"), function(r) r[2]))

pheno.dat <- read.csv('/brains/SWU4_phenotypic_data.csv')
pheno.dat$AGE_AT_SCAN_1 <- as.numeric(as.character(pheno.dat$AGE_AT_SCAN_1))
pheno.dat <- pheno.dat[!duplicated(pheno.dat$SUBID),]
pheno.dat <- pheno.dat[, c("SUBID", "AGE_AT_SCAN_1", "SEX")]

matched.idx <- sapply(subids, function(x) {
  match <- which(x == pheno.dat$SUBID)
  if (length(match) > 0) {
    return(match)
  } else {
    return(NULL)
  }
})
retain.idx <- which(sapply(matched.idx, is.null))
if (length(retain.idx) > 0) {
  matched.idx <- matched.idx[-retain.idx]
}

Y <- pheno.dat$SEX[matched.idx]

mask <- readNIfTI("/data/MNI152_T1_1mm_brain_mask.nii.gz")
ndim <- prod(dim(mask))
tic <- Sys.time()
res <- mclapply(1:length(fnames), function(i) {
  print(i)
  img <- readNIfTI(paste("/brains/brains/", fnames[i], sep=""))
  X <- apply(img, MARGIN=c(4), function(x) {x[mask > 0]})
  gc()
  return(c(X))
}, mc.cores=mc.cores/2)
toc <- Sys.time()
runtime.mtx <- toc - tic
gc()
D <- fm.as.matrix(do.call(rbind, res))
rm(res)

train.idx <- sample(1:length(Y), size=ceil(length(Y)/2), replace=FALSE)
test.idx <- (1:length(Y))[-train.idx]

D.train <- fm.get.rows(D, train.idx); D.test <- fm.get.rows(D, test.idx)
Y.train <- Y[train.idx]; Y.test <- Y[test.idx]

proj.algs <- list(PCA=flashx.pca, LOL=flashx.lol, LDA=flashx.lrlda, RP=flashx.rp, CCA=flashx.lrcca)
classifier.algs <- list(LDA=lda, RF=randomForest)

res <- lapply(names(proj.algs), function(proj.name) {
  proj.alg <- proj.algs[[proj.name]]
  tic <- Sys.time()
  res <- do.call(proj.alg, list(X=D.train, Y=Y.train, r=90))

  op <- (list(X.train=res$Xr, X.test=flashx.embed(D.test, res$A)))
  xv.res <- lapply(c(10, 30, 50, 70, 90), function(r) {
    X.train <- op$X.train[,1:r]; X.test <- op$X.test[,1:r]
    return(lapply(names(classifier.algs), function(class.name) {
      tryCatch({
        class.alg <- classifier.algs[[class.name]]
        trained.classifier <- do.call(class.alg, list(X.train, as.factor(Y.train)))
        Y.hat <- predict(trained.classifier, X.test)
        if (class.name == "LDA") {
          Y.hat <- Y.hat$class
        }
        acc <- mean(Y.hat == Y.test)
        return(data.frame(Dataset=dataset, Algorithm=proj.name, Classifier=class.name,
                          Accuracy=acc, Misclassification=1-acc, r=r))},
        error=function(e) {return(NULL)})
    }) %>%
      bind_rows())
  }) %>%
    bind_rows()
  toc <- Sys.time()
  saveRDS(list(xv=xv.res, time=toc - tic), sprintf('/brains/%s_flashlol_%s.rds', dataset, proj.name))
})
names(res) <- names(proj.algs)

saveRDS(list(result=res, Y.train=Y.train, Y.test=Y.test), sprintf('/brains/%s_flashlol.rds'))

