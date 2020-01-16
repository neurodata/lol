require(oro.nifti)
require(ramify)
require(parallel)
require(tidyverse)
require(profmem)
require(FlashR)

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

runtime.mtx <- toc - tic
gc()
D <- fm.as.matrix(do.call(rbind, res))
rm(res)

train.idx <- sample(1:length(Y), size=ceil(length(Y)/2), replace=FALSE)
test.idx <- (1:length(Y))[-train.idx]

D.train <- fm.get.rows(D, train.idx); D.test <- fm.get.rows(D, test.idx)
Y.train <- Y[train.idx]; Y.test <- Y[test.idx]

algs <- list(pca=flashx.pca, lol=flashx.lol, lrlda=flashx.lrlda, cca=flashx.lrcca, rp=flashx.rp)

res <- lapply(1:length(algs), function(i) {
  print(i)
  alg <- algs[[i]]
  tic <- Sys.time()
  res <- do.call(alg, list(X=D.train, Y=Y.train, r=91))
  toc <- Sys.time()
  res.in <- lapply(c(10, 30, 50, 70, 90), function(r) {
    tryCatch({
      return(list(X.train=flashx.embed(D.train, fm.get.cols(res$A, 1:r)),
                  X.test=flashx.embed(D.test, fm.get.cols(res$A, 1:r))))
    }, error=function(e) {return(NULL)})
  })
  return(list(xv=res.in, time=toc - tic))
})
names(res) <- names(algs)

saveRDS(list(result=res, Y.train=Y.train, Y.test=Y.test), '/brains/swu4_mini.rds')

