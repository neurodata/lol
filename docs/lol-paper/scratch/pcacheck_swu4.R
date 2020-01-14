require(oro.nifti)
require(ramify)
require(parallel)
require(tidyverse)
require(profmem)
require(flashR)

mc.cores <- detectCores()
fnames <- list.files('./')
d <- 100
mask <- readNIfTI("/data/MNI152_T1_1mm_brain_mask.nii.gz")
tic <- Sys.time()
D <- fm.as.matrix(do.call(rbind, mclapply(1:length(fnames), function(i) {
  print(i)
  img <- readNIfTI(fnames[i])
  X <- apply(img, MARGIN=c(4), function(x) {x[mask > 0]})
  X.fl <- c(X)
  return(X.fl)
}, mc.cores=mc.cores)))
toc <- Sys.time()
runtime.mtx <- toc - tic


tic <- Sys.time()
sv <- fm.svd(D)
toc <- Sys.time()
runtime.svd <- toc - tic

saveRDS(list(runtime.mtx=runtime.mtx, runtime.svd=runtime.svd, svd=sv,
             pca=fm.inner.prod(D, sv$u[,1:d])), file="./results_svd.rds")
