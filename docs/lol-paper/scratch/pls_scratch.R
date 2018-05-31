
# compute the cutoff for the particular trial to get an approximate elbow
# by computing the smallest r with an associated lhat within 5%
# of the global minimum lhat
compute_cutoff <- function(rs, lhats, t=0.05) {
  rs <- rs[complete.cases(lhats) & complete.cases(rs)]; lhats <- lhats[complete.cases(lhats) & complete.cases(rs)]
  sr.ix <- sort(rs, decreasing=FALSE, index.return=TRUE)$ix
  # compute minimum value
  min.lhat <- min(lhats)
  # compute minimum value + 5%
  lhat.thresh <- (1 + t)*min.lhat
  # find which indices are all below this
  lhat.below <- which(lhats <= lhat.thresh)
  rs.below <- rs[lhat.below]; lhats.below <- lhats[lhat.below]
  tmin.ix <- min(rs.below, index.return=TRUE)
  return(list(r=rs.below[tmin.ix], lhat=lhats.below[tmin.ix]))
}

opath <- './data/real_data/lda'
classifier.name = 'lda'
dataset.name = 'pmlb'

results <- readRDS(file.path(opath, paste(classifier.name, "_", dataset.name, "_results.rds", sep="")))
results <- results[complete.cases(results$lhat),]
nan.mean <- function(x) {mean(x, na.rm=TRUE)}
results.means <- aggregate(lhat ~ exp + alg + r + d + n + K, data = results, FUN = nan.mean)
random.results <- aggregate(lhat ~ exp + alg, data=subset(results, alg == "RandomGuess"), FUN=mean)

algs <-  c("LOL", "QOQ", "PLS", "PLSOL", "PLSOLK", "CCA", "LDA", "PCA", "RP")
acols <- c("#00FF00", "#00FF00", "#990000", "#990000", "#990000", "#AAAA55", "#000099", "#000099", "#000099")
linestyle <- c("solid", "dashed", "solid", "dashed", "dotted", "solid", "solid", "dashed", "dotted")
names(linestyle) <- algs
names(algs) <- acols
names(acols) <- algs
shapes <- c(21, 24, 21, 24, 23, 23, 21, 24, 23)
names(shapes) <- algs
exp_names <- unique(results$exp)


#nan.median <- function(x) median(x, na.rm=TRUE)
#results.medians <- aggregate(lhat ~ exp + alg + r + d + n + K, data = results, FUN = nan.median)
results.optimalr <- data.frame(exp=c(), alg=c(), r=c(), lhat=c(), fold=c())
plot.normlol.results <- data.frame(r=c(), lhat=c(), exp=c(), alg=c())
for (i in 1:length(exp_names)) {
  r.max <- max(results.means[results.means$exp == exp_names[i],]$r)
  lhat.mean <- random.results[random.results$exp == exp_names[i],]$lhat
  for (j in 1:length(algs)) {
    tryCatch({
      alg <- as.character(algs[j])
      ss <- results.means[results.means$exp == exp_names[i] & results.means$alg == algs[j],]
      rs <- ss$r; lhats <- ss$lhat
      min.result <- compute_cutoff(rs, lhats)
      r.min <- min.result$r; lhat.min <- min.result$lhat
      if (alg == 'LOL') {
        norm.r <- r.min
        norm.lhat <- lhat.min
      }
      #if (norm.r == 0) {
      #  if (r.min == 0) {
      #    r.rat <- 1
      #  } else {
      #    r.rat <- 10
      #  }
      #} else {
      #  r.rat <- (r.min - norm.r)/r.max
      #}
      r.rat <- (r.min - norm.r)/r.max
      #if (norm.lhat == 0) {
      #  if (lhat.min == 0) {
      #    lhat.rat <- 1
      #  } else {
      #    lhat.rat <- 10
      #  }
      #} else {
      #  lhat.rat <- lhat.min - norm.lhat
      #}
      lhat.rat <- (lhat.min - norm.lhat)/lhat.mean
      exp.fold <- results[results$exp == exp_names[i] & results$alg == algs[j] & results$r == r.min,]
      results.optimalr <- rbind(results.optimalr, data.frame(r=r.min, lhat=exp.fold$lhat,
                                                             fold=exp.fold$fold, alg=alg, exp=exp_names[i]))
      plot.normlol.results <- rbind(plot.normlol.results, data.frame(r=r.rat, lhat=lhat.rat,
                                                                     exp=exp_names[i], alg=alg))
    }, error=function(e) {NaN}, warning=function(w) {NaN})
  }
}

i=1
j=3

i.ss <- results.optimalr[results.optimalr$alg == algs[i],]
j.ss <- results.optimalr[results.optimalr$alg == algs[j],]
cmp <- merge(i.ss, j.ss, by=c("exp", "fold"), all=TRUE)

failed_pls_ss <- cmp[is.na(cmp$lhat.y),]
failed_pls_red <- aggregate(lhat.x ~ exp, data=failed_pls_ss, FUN=nan.mean)
failed_dsets <- failed_pls_red$exp

