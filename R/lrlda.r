lrlda = function(data, labels) {
  unique_classes = unique(labels)
  # Compute prior probabilities
  priors = sapply(unique_classes, function(x) sum(labels==x)/length(labels))
  counts = sapply(unique_classes, function(x) sum(labels==x))
  # Compute class means
  class_means = lapply(unique_classes, function(x) colMeans(data[labels==x,]))
  overall_mean = colMeans(data)

  sb = matrix(0, dim(data)[2], dim(data)[2])
  sw = matrix(0, dim(data)[2], dim(data)[2])

  for (i in 1:length(unique_classes)) {
    sb = sb + counts[i] * ((class_means[[i]] - overall_mean) %*% t(class_means[[i]] - overall_mean))
    x = data[labels==unique_classes[i],]
    sw = sw + t(x-class_means[[i]]) %*% (x-class_means[[i]])
  }

  w = eigen(ginv(sw) %*% sb)
}
