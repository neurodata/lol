require(oro.nifti)
require(ramify)
require(parallel)
require(tidyverse)
require(profmem)
require(FlashR)

mc.cores <- detectCores()
fnames <- list.files('/brains')
fnames <- fnames[fnames != "lost+found"]
d <- 100
mask <- readNIfTI("/data/MNI152_T1_1mm_brain_mask.nii.gz")
ndim <- prod(dim(mask))
tic <- Sys.time()
D <- array(NaN, dim=c(length(fnames, ndim)))
res <- mclapply(1:length(fnames), function(i) {
  print(i)
  img <- readNIfTI(paste("/brains/", fnames[i], sep=""))
  X <- apply(img, MARGIN=c(4), function(x) {x[mask > 0]})
  gc()
  return(c(X))
}, mc.cores=mc.cores/2)
toc <- Sys.time()
runtime.mtx <- toc - tic
gc()
D <- fm.as.matrix(do.call(rbind, res))
rm(res)
gc()

tic <- Sys.time()
sv <- fm.svd(D, nu=0, nv=d)
toc <- Sys.time()
pca=as.matrix(D %*% sv$v)
runtime.svd <- toc - tic

saveRDS(list(runtime.mtx=runtime.mtx, runtime.svd=runtime.svd, pca=pca), file="./test_svd_swu4.rds")
